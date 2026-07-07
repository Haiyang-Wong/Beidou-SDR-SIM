#include "bds_sdr.h"
#include <algorithm>
#include <atomic>
#include <cerrno>
#include <cstdlib>
#include <condition_variable>
#include <csignal>
#include <memory>
#include <mutex>
#include <thread>
#include <vector>

#ifdef BDS_SDR_WITH_UHD
#include <uhd/device.hpp>
#include <uhd/stream.hpp>
#include <uhd/usrp/multi_usrp.hpp>
#endif


namespace
{
constexpr int BDS_CODE_PHASE_BITS = 40;
constexpr uint64_t BDS_CODE_PHASE_SCALE = 1ULL << BDS_CODE_PHASE_BITS;
constexpr uint64_t BDS_CODE_PHASE_PERIOD = static_cast<uint64_t>(CA_SEQ_LEN_B1C) * BDS_CODE_PHASE_SCALE;
constexpr int BDS_CARR_PHASE_BITS = 48;
constexpr int64_t BDS_CARR_PHASE_SCALE = 1LL << BDS_CARR_PHASE_BITS;

uint64_t toCodePhaseFixed(double phase)
{
    return static_cast<uint64_t>(std::llround(phase * static_cast<double>(BDS_CODE_PHASE_SCALE)));
}

uint64_t toCodeStepFixed(double step)
{
    return static_cast<uint64_t>(std::llround(step * static_cast<double>(BDS_CODE_PHASE_SCALE)));
}

int64_t toCarrierPhaseFixed(double phase)
{
    return static_cast<int64_t>(std::llround(phase * static_cast<double>(BDS_CARR_PHASE_SCALE)));
}

int64_t toCarrierStepFixed(double step)
{
    return static_cast<int64_t>(std::llround(step * static_cast<double>(BDS_CARR_PHASE_SCALE)));
}

int64_t normalizeCarrierPhaseFixed(int64_t phase)
{
    phase %= BDS_CARR_PHASE_SCALE;
    if (phase < 0)
        phase += BDS_CARR_PHASE_SCALE;
    return phase;
}

int64_t wrapCarrierPhaseFixed(int64_t phase)
{
    if (phase >= BDS_CARR_PHASE_SCALE)
        return phase - BDS_CARR_PHASE_SCALE;
    if (phase < 0)
        return phase + BDS_CARR_PHASE_SCALE;
    return phase;
}

int carrierTableIndex(int64_t phase)
{
    return static_cast<int>((phase * 511) / BDS_CARR_PHASE_SCALE) & 511;
}

int codeIndex(uint64_t phase, int scale)
{
    return static_cast<int>((phase * static_cast<uint64_t>(scale)) >> BDS_CODE_PHASE_BITS);
}

#ifdef BDS_SDR_WITH_UHD
constexpr size_t BDS_SAMPLES_PER_TX_BUFFER = 64 * 1024;
constexpr double BDS_REALTIME_FIFO_MARGIN_SECONDS = 2.0;
constexpr const char *BDS_UHD_NUM_SEND_FRAMES = "1024";
volatile std::sig_atomic_t g_stop_requested = 0;

void sigintHandler(int)
{
    g_stop_requested = 1;
}

const char *txAsyncEventName(uhd::async_metadata_t::event_code_t event_code)
{
    switch (event_code)
    {
    case uhd::async_metadata_t::EVENT_CODE_BURST_ACK:
        return "BURST_ACK";
    case uhd::async_metadata_t::EVENT_CODE_UNDERFLOW:
        return "UNDERFLOW";
    case uhd::async_metadata_t::EVENT_CODE_SEQ_ERROR:
        return "SEQ_ERROR";
    case uhd::async_metadata_t::EVENT_CODE_TIME_ERROR:
        return "TIME_ERROR";
    case uhd::async_metadata_t::EVENT_CODE_UNDERFLOW_IN_PACKET:
        return "UNDERFLOW_IN_PACKET";
    case uhd::async_metadata_t::EVENT_CODE_SEQ_ERROR_IN_BURST:
        return "SEQ_ERROR_IN_BURST";
    case uhd::async_metadata_t::EVENT_CODE_USER_PAYLOAD:
        return "USER_PAYLOAD";
    }
    return "UNKNOWN";
}

class RealtimeIqFifo
{
public:
    explicit RealtimeIqFifo(size_t capacity_samples)
        : buffer_(capacity_samples * 2), capacity_samples_(capacity_samples)
    {
    }

    bool write(const short *samples, size_t sample_count)
    {
        std::unique_lock<std::mutex> lock(mutex_);
        write_ready_.wait(lock, [&]() {
            return stopped_ || finished_ || (capacity_samples_ - size_samples_) >= sample_count;
        });

        if (stopped_ || finished_)
            return false;

        copyIntoRing(samples, sample_count);
        size_samples_ += sample_count;
        lock.unlock();
        read_ready_.notify_one();
        return true;
    }

    size_t read(short *samples, size_t max_sample_count, bool *last_block = nullptr)
    {
        std::unique_lock<std::mutex> lock(mutex_);
        read_ready_.wait(lock, [&]() {
            return stopped_ || finished_ || size_samples_ > 0;
        });

        if (size_samples_ == 0)
        {
            if (last_block != nullptr)
                *last_block = true;
            return 0;
        }

        const size_t sample_count = std::min(max_sample_count, size_samples_);
        copyFromRing(samples, sample_count);
        size_samples_ -= sample_count;
        if (last_block != nullptr)
            *last_block = finished_ && size_samples_ == 0;
        lock.unlock();
        write_ready_.notify_one();
        return sample_count;
    }

    bool waitForSamples(size_t sample_count)
    {
        std::unique_lock<std::mutex> lock(mutex_);
        read_ready_.wait(lock, [&]() {
            return stopped_ || size_samples_ >= sample_count || (finished_ && size_samples_ > 0);
        });

        return size_samples_ > 0;
    }

    size_t size() const
    {
        std::lock_guard<std::mutex> lock(mutex_);
        return size_samples_;
    }

    void finish()
    {
        std::lock_guard<std::mutex> lock(mutex_);
        finished_ = true;
        read_ready_.notify_all();
        write_ready_.notify_all();
    }

    void stop()
    {
        std::lock_guard<std::mutex> lock(mutex_);
        stopped_ = true;
        read_ready_.notify_all();
        write_ready_.notify_all();
    }

private:
    void copyIntoRing(const short *samples, size_t sample_count)
    {
        const size_t samples_until_wrap = capacity_samples_ - head_samples_;
        const size_t first_samples = std::min(sample_count, samples_until_wrap);
        memcpy(&buffer_[head_samples_ * 2], samples, first_samples * 2 * sizeof(short));

        const size_t remaining_samples = sample_count - first_samples;
        if (remaining_samples > 0)
            memcpy(&buffer_[0], samples + first_samples * 2, remaining_samples * 2 * sizeof(short));

        head_samples_ = (head_samples_ + sample_count) % capacity_samples_;
    }

    void copyFromRing(short *samples, size_t sample_count)
    {
        const size_t samples_until_wrap = capacity_samples_ - tail_samples_;
        const size_t first_samples = std::min(sample_count, samples_until_wrap);
        memcpy(samples, &buffer_[tail_samples_ * 2], first_samples * 2 * sizeof(short));

        const size_t remaining_samples = sample_count - first_samples;
        if (remaining_samples > 0)
            memcpy(samples + first_samples * 2, &buffer_[0], remaining_samples * 2 * sizeof(short));

        tail_samples_ = (tail_samples_ + sample_count) % capacity_samples_;
    }

    std::vector<short> buffer_;
    size_t capacity_samples_;
    size_t head_samples_ = 0;
    size_t tail_samples_ = 0;
    size_t size_samples_ = 0;
    bool finished_ = false;
    bool stopped_ = false;
    mutable std::mutex mutex_;
    std::condition_variable read_ready_;
    std::condition_variable write_ready_;
};

struct RealtimeTxContext
{
    explicit RealtimeTxContext(size_t fifo_capacity_samples)
        : fifo(fifo_capacity_samples), buffer(BDS_SAMPLES_PER_TX_BUFFER * 2)
    {
    }

    RealtimeIqFifo fifo;
    uhd::usrp::multi_usrp::sptr usrp;
    uhd::tx_streamer::sptr stream;
    uhd::tx_metadata_t metadata;
    std::vector<short> buffer;
    std::thread thread;
    size_t prebuffer_samples = 0;
    size_t fifo_blocks = 0;
    size_t prebuffer_blocks = 0;
    std::atomic<long long> samples_sent{0};
    std::mutex error_mutex;
    std::string error;
};

void setRealtimeError(RealtimeTxContext *ctx, const std::string &message)
{
    std::lock_guard<std::mutex> lock(ctx->error_mutex);
    if (ctx->error.empty())
        ctx->error = message;
}

std::string getRealtimeError(RealtimeTxContext *ctx)
{
    std::lock_guard<std::mutex> lock(ctx->error_mutex);
    return ctx->error;
}

void realtimeTxTask(RealtimeTxContext *ctx)
{
    try
    {
        if (!ctx->fifo.waitForSamples(ctx->prebuffer_samples))
            return;

        cout << "USRP TX      : prebuffered " << ctx->fifo.size() << " samples, starting stream" << endl;

        bool first_packet = true;
        bool end_of_burst_sent = false;
        while (!g_stop_requested)
        {
            bool last_packet = false;
            const size_t sample_count = ctx->fifo.read(ctx->buffer.data(), BDS_SAMPLES_PER_TX_BUFFER, &last_packet);
            if (sample_count == 0)
                break;

            uhd::tx_metadata_t metadata = ctx->metadata;
            metadata.start_of_burst = first_packet;
            metadata.end_of_burst = last_packet;
            first_packet = false;
            end_of_burst_sent = last_packet;

            const size_t sent = ctx->stream->send(ctx->buffer.data(), sample_count, metadata, 1.0);
            ctx->samples_sent += static_cast<long long>(sent);
            if (sent != sample_count)
            {
                std::ostringstream oss;
                oss << "USRP sent " << sent << " samples, expected " << sample_count;
                setRealtimeError(ctx, oss.str());
                ctx->fifo.stop();
                break;
            }
        }

        if (end_of_burst_sent)
        {
            uhd::async_metadata_t async_metadata{};
            bool burst_ack_received = false;
            for (int attempt = 0; attempt < 50; ++attempt)
            {
                if (!ctx->stream->recv_async_msg(async_metadata, 0.1))
                    continue;

                if (async_metadata.event_code == uhd::async_metadata_t::EVENT_CODE_BURST_ACK)
                {
                    burst_ack_received = true;
                    break;
                }

                cout << "\nUSRP TX async: " << txAsyncEventName(async_metadata.event_code) << endl;
            }

            if (!burst_ack_received)
                cout << "\nUSRP TX async: no BURST_ACK received before shutdown" << endl;
        }
    }
    catch (const std::exception &ex)
    {
        setRealtimeError(ctx, ex.what());
        ctx->fifo.stop();
    }
}

std::unique_ptr<RealtimeTxContext> startRealtimeTx(const option_t &opt, size_t generated_block_samples)
{
    g_stop_requested = 0;
    signal(SIGINT, sigintHandler);

    const double block_seconds = static_cast<double>(generated_block_samples) / opt.samp_rate;
    const size_t prebuffer_blocks = std::max<size_t>(
        1,
        static_cast<size_t>(std::ceil(opt.tx_prebuffer / block_seconds)));
    const size_t margin_blocks = std::max<size_t>(
        1,
        static_cast<size_t>(std::ceil(BDS_REALTIME_FIFO_MARGIN_SECONDS / block_seconds)));
    const size_t fifo_blocks = prebuffer_blocks + margin_blocks;

    auto ctx = std::make_unique<RealtimeTxContext>(generated_block_samples * fifo_blocks);
    ctx->prebuffer_blocks = prebuffer_blocks;
    ctx->fifo_blocks = fifo_blocks;
    ctx->prebuffer_samples = generated_block_samples * prebuffer_blocks;

    uhd::device_addr_t device_addr(std::string(opt.usrp_args));
    if (std::string(opt.usrp_args).empty())
    {
        cout << "USRP device  : auto-detecting first available UHD device" << endl;
        const uhd::device_addrs_t devices = uhd::device::find(uhd::device_addr_t());
        if (devices.empty())
            throw std::runtime_error("No UHD devices found. Check USB/network connection and uhd_find_devices.");

        device_addr = devices.front();
    }

    if (!device_addr.has_key("num_send_frames"))
        device_addr["num_send_frames"] = BDS_UHD_NUM_SEND_FRAMES;

    cout << "USRP device  : " << device_addr.to_string() << endl;

    ctx->usrp = uhd::usrp::multi_usrp::make(device_addr);
    ctx->usrp->set_clock_source("internal");
    ctx->usrp->set_time_source("internal");

    ctx->usrp->set_tx_rate(opt.samp_rate);
    ctx->usrp->set_tx_freq(opt.tx_freq);
    ctx->usrp->set_tx_gain(opt.tx_gain);

    cout << "Actual TX rate : " << ctx->usrp->get_tx_rate() / 1e6 << " MHz" << endl;
    cout << "Actual TX freq : " << ctx->usrp->get_tx_freq() / 1e6 << " MHz" << endl;
    cout << "Actual TX gain : " << ctx->usrp->get_tx_gain() << " dB" << endl;
    cout << "TX antenna     : " << ctx->usrp->get_tx_antenna() << endl;

    uhd::stream_args_t stream_args("sc16", "sc16");
    ctx->stream = ctx->usrp->get_tx_stream(stream_args);

    ctx->metadata.start_of_burst = false;
    ctx->metadata.end_of_burst = false;
    ctx->metadata.has_time_spec = false;

    ctx->thread = std::thread(realtimeTxTask, ctx.get());
    cout << "USRP TX      : FIFO capacity=" << ctx->fifo_blocks
         << " blocks, prebuffer=" << ctx->prebuffer_blocks
         << " blocks (" << fixed << setprecision(1) << ctx->prebuffer_blocks * block_seconds
         << " s)" << defaultfloat << setprecision(6) << endl;
    cout << "USRP realtime TX task created. Press Ctrl+C to stop streaming." << endl;
    return ctx;
}
#endif

struct TrajectoryPoint
{
    double time;
    double llh[3];
};

std::string trimCopy(const std::string &value)
{
    const size_t first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos)
        return "";

    const size_t last = value.find_last_not_of(" \t\r\n");
    return value.substr(first, last - first + 1);
}

bool parseDoubleStrict(const std::string &token, double &value)
{
    const std::string trimmed = trimCopy(token);
    if (trimmed.empty())
        return false;

    errno = 0;
    char *end = nullptr;
    value = std::strtod(trimmed.c_str(), &end);
    return errno != ERANGE && end == trimmed.c_str() + trimmed.size() && std::isfinite(value);
}

std::vector<std::string> splitCsvLine(const std::string &line)
{
    std::vector<std::string> fields;
    std::stringstream ss(line);
    std::string field;

    while (std::getline(ss, field, ','))
        fields.push_back(field);

    if (!line.empty() && line.back() == ',')
        fields.emplace_back();

    return fields;
}

bool parseTrajectoryRow(const std::string &line, TrajectoryPoint &point, std::string &error)
{
    const std::vector<std::string> fields = splitCsvLine(line);
    if (fields.size() != 4)
    {
        error = "expected 4 CSV fields: time,lat,lon,alt";
        return false;
    }

    if (!parseDoubleStrict(fields[0], point.time) ||
        !parseDoubleStrict(fields[1], point.llh[0]) ||
        !parseDoubleStrict(fields[2], point.llh[1]) ||
        !parseDoubleStrict(fields[3], point.llh[2]))
    {
        error = "expected numeric time,lat,lon,alt fields";
        return false;
    }

    if (point.time < 0.0)
    {
        error = "time must be non-negative";
        return false;
    }
    if (point.llh[0] < -90.0 || point.llh[0] > 90.0)
    {
        error = "latitude must be in [-90,90] degrees";
        return false;
    }
    if (point.llh[1] < -180.0 || point.llh[1] > 360.0)
    {
        error = "longitude must be in [-180,360] degrees";
        return false;
    }

    return true;
}

bool loadTrajectoryCsv(const char *path, std::vector<TrajectoryPoint> &points)
{
    std::ifstream input(path);
    if (!input.is_open())
    {
        cerr << "ERROR: Failed to open trajectory CSV: " << path << endl;
        return false;
    }

    std::string line;
    int line_number = 0;
    bool skipped_header = false;

    while (std::getline(input, line))
    {
        line_number++;
        if (line_number == 1 && line.size() >= 3 &&
            static_cast<unsigned char>(line[0]) == 0xEF &&
            static_cast<unsigned char>(line[1]) == 0xBB &&
            static_cast<unsigned char>(line[2]) == 0xBF)
        {
            line.erase(0, 3);
        }

        line = trimCopy(line);
        if (line.empty())
            continue;

        TrajectoryPoint point{};
        std::string error;
        if (!parseTrajectoryRow(line, point, error))
        {
            if (points.empty() && !skipped_header)
            {
                skipped_header = true;
                continue;
            }

            cerr << "ERROR: Invalid trajectory CSV line " << line_number << ": " << error << endl;
            return false;
        }

        if (points.empty())
        {
            if (std::fabs(point.time) > 1e-3)
            {
                cerr << "ERROR: First trajectory time must be 0 or close to 0 seconds." << endl;
                return false;
            }
        }
        else if (point.time <= points.back().time)
        {
            cerr << "ERROR: Trajectory time must be strictly increasing at line " << line_number << "." << endl;
            return false;
        }

        points.push_back(point);
    }

    if (points.empty())
    {
        cerr << "ERROR: Trajectory CSV has no valid trajectory points: " << path << endl;
        return false;
    }

    return true;
}

void llhDegreesToXyz(const double llh_degrees[3], double xyz[3])
{
    double llh_radians[3] = {
        llh_degrees[0] / R2D,
        llh_degrees[1] / R2D,
        llh_degrees[2]};
    llh2xyz(llh_radians, xyz);
}

void interpolateTrajectoryXyz(
    const std::vector<TrajectoryPoint> &points,
    double time,
    size_t &segment_index,
    double xyz[3])
{
    double llh_degrees[3] = {0.0, 0.0, 0.0};

    if (time <= points.front().time)
    {
        std::copy(points.front().llh, points.front().llh + 3, llh_degrees);
    }
    else if (time >= points.back().time)
    {
        std::copy(points.back().llh, points.back().llh + 3, llh_degrees);
    }
    else
    {
        while (segment_index + 1 < points.size() && points[segment_index + 1].time < time)
            segment_index++;

        const TrajectoryPoint &start = points[segment_index];
        const TrajectoryPoint &end = points[segment_index + 1];
        const double ratio = (time - start.time) / (end.time - start.time);
        for (int i = 0; i < 3; ++i)
            llh_degrees[i] = start.llh[i] + ratio * (end.llh[i] - start.llh[i]);
    }

    llhDegreesToXyz(llh_degrees, xyz);
}

void receiverXyzAtTime(
    const option_t &opt,
    const std::vector<TrajectoryPoint> &trajectory,
    double time,
    size_t &trajectory_segment_index,
    double xyz[3])
{
    if (opt.use_trajectory)
        interpolateTrajectoryXyz(trajectory, time, trajectory_segment_index, xyz);
    else
        llhDegreesToXyz(opt.llh, xyz);
}

}

void *runBdsTask(sim_t *sim)
{

    cout << "BDS SDR Simulation Task Started." << endl;
    sim->status = 0;

    int duration;                             // 仿真持续时间，单位为1秒
    vector<ephem_t> eph_vector[MAX_SAT];       // 存储所有卫星的星历数据
    double xyz[3];
    std::vector<TrajectoryPoint> trajectory;
    size_t trajectory_segment_index = 0;
    int sv;                                    // 卫星prn编号
    ephem_t eph;                               // 临时存储星历数据
    bdstime_t bdt_min, bdt_max;                // 临时存储BDT时间
    datetime_t t0, tmin, tmax;                 // 临时存储日期时间数据
    bdstime_t g0;                              // 仿真起始时间
    bdstime_t brx;                             // 接收机接收信号的时间
    vector<int> current_eph_index(MAX_SAT, 0); // 当前使用的星历索引
    double dt = 0.1;           // 时间增量，单位为秒
    channel_t channels[MAX_CHAN];              // 信道数组
    double elvmask = sim->opt.elv_mask;        // 仰角掩模，单位为度
    vector<int> allocated_channels(MAX_SAT);   // allocated channel for each prn
    ephem_t bds_almanac[63] = {};          // B1C almanac used by subframe 3
    int iq_buff_size;                          // IQ缓冲区大小
    double delt_between_samples;               // 采样点间的相位增量
    short *iq_buff = NULL;                     // I/Q缓冲区指针
    FILE *fp = NULL;                           // 标准C文件指针，用于文件输入输出操作
    char out_file[MAX_CHAR];
    int num_itime; // 当前仿真时间步数

    duration = sim->opt.duration;
    if (sim->opt.use_trajectory)
    {
        if (!loadTrajectoryCsv(sim->opt.traj_file, trajectory))
        {
            sim->status = 1;
            return nullptr;
        }
        cout << "Trajectory   : loaded " << trajectory.size() << " points" << endl;
    }

    receiverXyzAtTime(sim->opt, trajectory, 0.0, trajectory_segment_index, xyz);
    int eph_count = readBdsB1CEphemerisCpp(eph_vector, sim->opt.nav_file);
    computeBdtMinMax(eph_vector, &bdt_min, &bdt_max, &tmin, &tmax);
    setScenarioStartTimeBds(&sim->opt.g0, bdt_min, bdt_max, &t0, &tmin, &tmax, sim->opt.time_overwrite, eph_count, eph_vector);
    g0 = sim->opt.g0;

    brx = g0; // 接收机接收信号的时间从仿真起始时间开始
    initChannels(channels, allocated_channels);
    auto release_channels = [&channels]() {
        for (int i = 0; i < MAX_CHAN; i++)
        {
            free(channels[i].ca_b1c_data);
            free(channels[i].ca_b1c_pilot11);
            free(channels[i].ca_b1c_pilot61);
            free(channels[i].weighted_b1c_data);
            free(channels[i].weighted_b1c_pilot11);
            free(channels[i].weighted_b1c_pilot61);
            free(channels[i].b1c_pilot_sub_code);
            free(channels[i].nav_bit);
        }
    };

    auto refresh_b1c_almanac = [&]() {
        const ephem_t *eph_list[63] = {};

        for (int alm_sv = 0; alm_sv < 63; ++alm_sv)
        {
            if (eph_vector[alm_sv].empty())
                continue;

            int idx = matchEpoch(brx, eph_vector[alm_sv]);
            if (idx < 0)
                continue;

            current_eph_index[alm_sv] = idx;
            eph_list[alm_sv] = &eph_vector[alm_sv][idx];
        }

        memset(bds_almanac, 0, sizeof(bds_almanac));
        const int tow_seconds = static_cast<int>((static_cast<double>(brx.milliseconds) + brx.sub_milliseconds) /
                                                BDS_MILLISECONDS_IN_SECOND);
        completeBdsAlmanacFromEphem(eph_list, bds_almanac, brx.week, tow_seconds);
    };
    refresh_b1c_almanac();
    int nsat = allocateChannel(channels, eph_vector, brx, xyz, elvmask, current_eph_index, allocated_channels, bds_almanac);
    if (nsat <= 0)
    {
        cerr << "ERROR: No visible BDS satellites at the scenario start time." << endl;
        cerr << "       Check --utc, --lla, and --elv, or lower the elevation mask." << endl;
        release_channels();
        return nullptr;
    }

    iq_buff_size = (int)(sim->opt.samp_rate / 10);      // 就是0.1s内生成的采样点的数量
    delt_between_samples = 1.0 / sim->opt.samp_rate;    // 采样点间的时间间隔，单位为秒
    iq_buff = (short *)calloc(2 * iq_buff_size, sizeof(short));
    if (iq_buff == NULL)
    {
        cerr << "ERROR: Failed to allocate I/Q sample buffer." << endl;
        sim->status = 1;
        release_channels();
        return nullptr;
    }

    const unsigned int hw_threads = std::thread::hardware_concurrency();
    const int requested_workers = static_cast<int>(hw_threads == 0 ? 4 : hw_threads);
    const int render_worker_count = std::max(1, std::min({requested_workers, 8, iq_buff_size}));

#ifdef BDS_SDR_WITH_UHD
    std::unique_ptr<RealtimeTxContext> realtime_tx;
    if (sim->opt.use_usrp)
    {
        try
        {
            realtime_tx = startRealtimeTx(sim->opt, static_cast<size_t>(iq_buff_size));
        }
        catch (const std::exception &ex)
        {
            cerr << "ERROR: Failed to initialize USRP realtime TX: " << ex.what() << endl;
            sim->status = 1;
            free(iq_buff);
            release_channels();
            return nullptr;
        }
    }
#else
    if (sim->opt.use_usrp)
    {
        cerr << "ERROR: This build was configured without UHD, so --usrp is unavailable." << endl;
        cerr << "       Install UHD development libraries and reconfigure CMake to enable realtime TX." << endl;
        sim->status = 1;
        free(iq_buff);
        release_channels();
        return nullptr;
    }
#endif

    if (!sim->opt.use_usrp)
    {
        strcpy(out_file, sim->opt.out_file);
        if (strcmp("-", out_file))
        {
            if (NULL == (fp = fopen(out_file, "wb")))
            {
                fprintf(stderr, "ERROR: Failed to open output file.\n");
                exit(1);
            }
        }
        else
        {
            fp = stdout;
        }
    }
    cout << "------------------ Start generating simulation signal -------------------" << endl;
    num_itime = duration * 10;
    brx = addBdsTime(brx, dt);
    for (int itime = 1; itime < num_itime; itime++)
    {
#ifdef BDS_SDR_WITH_UHD
        if (sim->opt.use_usrp && g_stop_requested)
            break;
#endif
        std::ostringstream progress;
        progress << "Simulation time: " << fixed << setprecision(1) << itime * 0.1 << " s";
        cout << '\r' << progress.str() << flush;
        receiverXyzAtTime(sim->opt, trajectory, itime * dt, trajectory_segment_index, xyz);
        for (int i = 0; i < MAX_CHAN; i++)
        {
            if (channels[i].prn > 0)
            {
                range_t rho;
                sv = channels[i].prn - 1;
                eph = eph_vector [sv][current_eph_index[sv]];
                computeRange(&rho, eph, brx, xyz, channels[i].prn);
                channels[i].azel[0] = rho.azel[0]; // 方位角
                channels[i].azel[1] = rho.azel[1]; // 高度角
                computeCodePhase(&channels[i], rho, dt, brx);   
            }
        }
        struct ActiveChannelRenderState
        {
            int channel_index;
            const double *weighted_b1c_data;
            const double *weighted_b1c_pilot11;
            const double *weighted_b1c_pilot61;
            const short *b1c_pilot_sub_code;
            const short *nav_bit;
            uint64_t code_phase;
            uint64_t code_step;
            int64_t carr_phase;
            int64_t carr_step;
            int ibit;
        };

        std::vector<ActiveChannelRenderState> active_channels;
        active_channels.reserve(MAX_CHAN);
        auto buildActiveChannels = [&]() {
            std::vector<ActiveChannelRenderState> active;
            active.reserve(MAX_CHAN);
            for (int i_chan = 0; i_chan < MAX_CHAN; ++i_chan)
            {
                if (channels[i_chan].prn <= 0)
                    continue;

                ActiveChannelRenderState state{};
                state.channel_index = i_chan;
                state.weighted_b1c_data = channels[i_chan].weighted_b1c_data;
                state.weighted_b1c_pilot11 = channels[i_chan].weighted_b1c_pilot11;
                state.weighted_b1c_pilot61 = channels[i_chan].weighted_b1c_pilot61;
                state.b1c_pilot_sub_code = channels[i_chan].b1c_pilot_sub_code;
                state.nav_bit = channels[i_chan].nav_bit;
                state.code_phase = toCodePhaseFixed(channels[i_chan].code_phase);
                state.code_step = toCodeStepFixed(channels[i_chan].f_code * delt_between_samples);
                state.carr_phase = toCarrierPhaseFixed(channels[i_chan].carr_phase);
                state.carr_step = toCarrierStepFixed(channels[i_chan].f_carr * delt_between_samples);
                state.ibit = channels[i_chan].ibit;
                active.push_back(state);
            }
            return active;
        };

        auto renderSegment = [&](const std::vector<ActiveChannelRenderState> &segment_channels,
                                 int segment_begin,
                                 int segment_end) {
            if (segment_begin >= segment_end)
                return;
            if (segment_channels.empty())
            {
                memset(iq_buff + segment_begin * 2, 0, static_cast<size_t>(segment_end - segment_begin) * 2 * sizeof(short));
                return;
            }

            const int segment_samples = segment_end - segment_begin;
            const int worker_count = std::min(render_worker_count, segment_samples);
            const int chunk_size = (segment_samples + worker_count - 1) / worker_count;
            std::vector<std::thread> workers;
            workers.reserve(worker_count);

            for (int worker = 0; worker < worker_count; ++worker)
            {
                const int sample_begin = segment_begin + worker * chunk_size;
                const int sample_end = std::min(segment_end, sample_begin + chunk_size);
                if (sample_begin >= sample_end)
                    break;

                workers.emplace_back([&, sample_begin, sample_end, segment_begin]() {
                    std::vector<ActiveChannelRenderState> local_channels = segment_channels;
                    const int segment_offset = sample_begin - segment_begin;
                    for (auto &state : local_channels)
                    {
                        uint64_t code_phase = state.code_phase + static_cast<uint64_t>(segment_offset) * state.code_step;
                        int code_wraps = static_cast<int>(code_phase / BDS_CODE_PHASE_PERIOD);
                        state.code_phase = code_phase - static_cast<uint64_t>(code_wraps) * BDS_CODE_PHASE_PERIOD;
                        state.ibit += code_wraps;

                        int64_t carr_phase = state.carr_phase + static_cast<int64_t>(segment_offset) * state.carr_step;
                        state.carr_phase = normalizeCarrierPhaseFixed(carr_phase);
                    }

                    for (int isample = sample_begin; isample < sample_end; ++isample)
                    {
                        int i_qua = 0;
                        int q_qua = 0;

                        for (auto &state : local_channels)
                        {
                            if (state.code_phase >= BDS_CODE_PHASE_PERIOD)
                            {
                                state.code_phase -= BDS_CODE_PHASE_PERIOD;
                                state.ibit++;
                            }

                            const int table_index = carrierTableIndex(state.carr_phase);
                            const int cos_ph = cos_table_512[table_index];
                            const int sin_ph = sin_table_512[table_index];
                            const int code_index_boc11 = codeIndex(state.code_phase, 2);
                            const int code_index_boc61 = codeIndex(state.code_phase, 12);

                            const double data_chip = state.weighted_b1c_data[code_index_boc11];
                            const double pilot11_chip = state.weighted_b1c_pilot11[code_index_boc11];
                            const double pilot61_chip = state.weighted_b1c_pilot61[code_index_boc61];
                            const int data_bit = state.nav_bit[state.ibit];
                            const int secondary_code = state.b1c_pilot_sub_code[state.ibit];
                            const double A = data_chip * data_bit + pilot61_chip * secondary_code;
                            const double B = pilot11_chip * secondary_code;

                            i_qua += A * cos_ph - B * sin_ph;
                            q_qua += A * sin_ph + B * cos_ph;
                            state.code_phase += state.code_step;
                            state.carr_phase = wrapCarrierPhaseFixed(state.carr_phase + state.carr_step);
                        }

                        iq_buff[isample * 2] = static_cast<short>(i_qua);
                        iq_buff[isample * 2 + 1] = static_cast<short>(q_qua);
                    }
                });
            }

            for (std::thread &worker : workers)
                worker.join();
        };

        auto advanceActiveChannels = [&](std::vector<ActiveChannelRenderState> &segment_channels, int sample_count) {
            if (sample_count <= 0)
                return;

            const int last_sample = sample_count - 1;
            for (auto &state : segment_channels)
            {
                const uint64_t last_code_phase = state.code_phase + static_cast<uint64_t>(last_sample) * state.code_step;
                const int applied_code_wraps = static_cast<int>(last_code_phase / BDS_CODE_PHASE_PERIOD);
                const uint64_t final_code_phase = state.code_phase + static_cast<uint64_t>(sample_count) * state.code_step -
                                                  static_cast<uint64_t>(applied_code_wraps) * BDS_CODE_PHASE_PERIOD;
                state.code_phase = final_code_phase;
                state.ibit += applied_code_wraps;
                state.carr_phase = normalizeCarrierPhaseFixed(state.carr_phase + static_cast<int64_t>(sample_count) * state.carr_step);
            }
        };

        auto nextNavBoundaryOffset = [&](const std::vector<ActiveChannelRenderState> &segment_channels,
                                         int samples_remaining) {
            int next_boundary = samples_remaining;
            if (samples_remaining <= 0)
                return next_boundary;

            const int last_sample = samples_remaining - 1;
            for (const auto &state : segment_channels)
            {
                if (state.code_step == 0 || state.ibit >= NAV_SEQ_LEN_B1C)
                    continue;

                const uint64_t max_code_phase = state.code_phase + static_cast<uint64_t>(last_sample) * state.code_step;
                const int max_code_wraps = static_cast<int>(max_code_phase / BDS_CODE_PHASE_PERIOD);
                if (state.ibit + max_code_wraps < NAV_SEQ_LEN_B1C)
                    continue;

                const int wraps_to_nav = NAV_SEQ_LEN_B1C - state.ibit;
                const uint64_t target_phase = static_cast<uint64_t>(wraps_to_nav) * BDS_CODE_PHASE_PERIOD;
                int boundary = 0;
                if (state.code_phase < target_phase)
                {
                    const uint64_t delta = target_phase - state.code_phase;
                    boundary = static_cast<int>((delta + state.code_step - 1) / state.code_step);
                }
                next_boundary = std::min(next_boundary, boundary);
            }
            return next_boundary;
        };

        auto applyNavBoundaryAtCursor = [&](std::vector<ActiveChannelRenderState> &segment_channels, int cursor) {
            for (auto &state : segment_channels)
            {
                if (state.code_phase < BDS_CODE_PHASE_PERIOD || state.ibit + 1 < NAV_SEQ_LEN_B1C)
                    continue;

                state.code_phase -= BDS_CODE_PHASE_PERIOD;
                state.ibit = 0;

                channel_t &chan = channels[state.channel_index];
                chan.code_phase = static_cast<double>(state.code_phase) / static_cast<double>(BDS_CODE_PHASE_SCALE);
                chan.carr_phase = static_cast<double>(state.carr_phase) / static_cast<double>(BDS_CARR_PHASE_SCALE);
                chan.ibit = 0;

                sv = chan.prn - 1;
                eph = eph_vector[sv][current_eph_index[sv]];
                bdstime_t sample_rx_time = addBdsTime(brx, cursor * delt_between_samples);
                bdstime_t g_tx = computeSatelliteTxTime(sample_rx_time, chan.rho0.range);
                refresh_b1c_almanac();
                generateB1CNavMessage(g_tx, &chan, &eph, bds_almanac);
            }
        };

        auto storeActiveChannels = [&]() {
            for (const auto &state : active_channels)
            {
                channel_t &chan = channels[state.channel_index];
                chan.code_phase = static_cast<double>(state.code_phase) / static_cast<double>(BDS_CODE_PHASE_SCALE);
                chan.carr_phase = static_cast<double>(state.carr_phase) / static_cast<double>(BDS_CARR_PHASE_SCALE);
                chan.ibit = state.ibit;
            }
        };

        active_channels = buildActiveChannels();
        int segment_begin = 0;
        while (segment_begin < iq_buff_size)
        {
            const int boundary_offset = nextNavBoundaryOffset(active_channels, iq_buff_size - segment_begin);
            const int segment_end = segment_begin + boundary_offset;
            renderSegment(active_channels, segment_begin, segment_end);
            advanceActiveChannels(active_channels, boundary_offset);
            segment_begin = segment_end;

            if (segment_begin < iq_buff_size)
                applyNavBoundaryAtCursor(active_channels, segment_begin);
        }
        storeActiveChannels();

        if (!sim->opt.use_usrp)
            fwrite(iq_buff, sizeof(short), 2 * iq_buff_size, fp);

#ifdef BDS_SDR_WITH_UHD
        if (realtime_tx)
        {
            const std::string tx_error = getRealtimeError(realtime_tx.get());
            if (!tx_error.empty())
            {
                cerr << "\nERROR: USRP realtime TX failed: " << tx_error << endl;
                sim->status = 1;
                break;
            }

            if (!realtime_tx->fifo.write(iq_buff, static_cast<size_t>(iq_buff_size)))
                break;
        }
#endif
        int ibrx = (int)(((static_cast<double>(brx.milliseconds) + brx.sub_milliseconds) / 100.0) + 0.5);
        if ((int)fmodf(ibrx, 300) == 0)
        {
            refresh_b1c_almanac();
            allocateChannel(channels, eph_vector, brx, xyz, elvmask, current_eph_index, allocated_channels, bds_almanac);
        }
        brx = addBdsTime(brx, dt);

    }
#ifdef BDS_SDR_WITH_UHD
    if (realtime_tx)
    {
        realtime_tx->fifo.finish();
        if (realtime_tx->thread.joinable())
            realtime_tx->thread.join();

        const std::string tx_error = getRealtimeError(realtime_tx.get());
        if (!tx_error.empty())
        {
            cerr << "\nERROR: USRP realtime TX failed: " << tx_error << endl;
            sim->status = 1;
        }
        cout << "\nTotal USRP samples sent: " << realtime_tx->samples_sent.load() << endl;
    }
#endif
    cout << endl;
    cout << "------------------ Simulation finished -------------------" << endl;
    if (fp != NULL && fp != stdout)
        fclose(fp);
    free(iq_buff);
    release_channels();
    return nullptr;
}

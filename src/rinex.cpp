#include "bds_sdr.h"
static std::string convertD2E(std::string line)
{
    std::replace(line.begin(), line.end(), 'D', 'E');
    std::replace(line.begin(), line.end(), 'd', 'E');
    return line;
}

static double readFixedDoubleField(const std::string &line, size_t offset, size_t width = 19)
{
    if (offset >= line.size())
        return 0.0;

    const std::string field = line.substr(offset, std::min(width, line.size() - offset));
    if (field.find_first_not_of(' ') == std::string::npos)
        return 0.0;

    std::istringstream iss(field);
    double value = 0.0;
    iss >> value;
    return iss.fail() ? 0.0 : value;
}

int readContentsData(const std::string &line, double data[39], int data_offset, datetime_t *time)
{
    const std::string normalized_line = convertD2E(line);
    int svid = 0;

    if (time != nullptr)
    {
        std::istringstream time_stream(normalized_line.substr(4, 19));
        time_stream >> time->y >> time->m >> time->d >> time->hh >> time->mm >> time->sec;

        if (normalized_line.size() > 2)
        {
            std::istringstream svid_stream(normalized_line.substr(1, 2));
            svid_stream >> svid;
        }

        data[data_offset] = readFixedDoubleField(normalized_line, 23);
        data[data_offset + 1] = readFixedDoubleField(normalized_line, 42);
        data[data_offset + 2] = readFixedDoubleField(normalized_line, 61);
        return svid;
    }

    data[data_offset] = readFixedDoubleField(normalized_line, 4);
    data[data_offset + 1] = readFixedDoubleField(normalized_line, 23);
    data[data_offset + 2] = readFixedDoubleField(normalized_line, 42);
    data[data_offset + 3] = readFixedDoubleField(normalized_line, 61);
    return svid;
}


double normalizeBdsAngle(double angle)
{
    while (angle < -K_PI)
        angle += 2.0 * K_PI;
    while (angle >= K_PI)
        angle -= 2.0 * K_PI;
    return angle;
}

int alignBdsAlmanacToa4096(int toa, int *week)
{
    toa = (toa + 2048) / 4096 * 4096;

    if (toa >= static_cast<int>(SECONDS_IN_WEEK))
    {
        toa -= static_cast<int>(SECONDS_IN_WEEK);
        if (week != nullptr)
            ++(*week);
    }

    return toa;
}

unsigned char getBdsSatTypeFromEphem(const ephem_t *eph)
{
    if (eph == nullptr)
        return 0;

    const double axis = (eph->axis > 0.0) ? eph->axis : eph->sqrt_a * eph->sqrt_a;

    if (eph->svid <= 5 || eph->svid >= 59)
        return 1;

    return (axis > 4.0e7) ? 2 : 3;
}

ephem_t deriveBdsAlmanacFromEphem(const ephem_t *eph, int alm_week, int alm_toa)
{
    ephem_t alm{};
    if (eph == nullptr)
        return alm;

    alm.valid = eph->valid & 1;
    alm.health = static_cast<unsigned char>(eph->health);
    alm.svid = eph->svid;

    if ((alm.valid & 1) == 0)
        return alm;

    alm.flag = getBdsSatTypeFromEphem(eph);
    alm.toa = alm_toa;
    alm.week = alm_week;

    const int dt = (alm_week - eph->week) * static_cast<int>(SECONDS_IN_WEEK) +
                   (alm_toa - eph->toe);
    const double axis = (eph->axis > 0.0) ? eph->axis : eph->sqrt_a * eph->sqrt_a;
    const double mean_motion = (eph->n != 0.0 || eph->sqrt_a <= 0.0 || axis <= 0.0)
                                   ? eph->n
                                   : CGCS2000_SQRT_GM / (eph->sqrt_a * axis) + eph->delta_n;

    alm.ecc = eph->ecc;
    alm.sqrt_a = eph->sqrt_a;
    alm.w = eph->w;
    alm.omega_dot = eph->omega_dot;
    alm.af1 = eph->af1;

    alm.m0 = normalizeBdsAngle(eph->m0 + mean_motion * dt);
    alm.omega0 = normalizeBdsAngle(eph->omega0 + eph->omega_dot * dt);
    alm.i0 = eph->i0 + eph->idot * dt;
    alm.af0 = eph->af0 + eph->af1 * dt;

    return alm;
}

void completeBdsAlmanacFromEphem(
    const ephem_t *eph_list[63],
    ephem_t alm_out[63],
    int current_bds_week,
    int current_bds_tow_seconds)
{
    int alm_toa = -1;
    int alm_week = current_bds_week;

    for (int i = 0; i < 63; ++i)
    {
        if ((alm_out[i].valid & 1) != 0)
        {
            alm_toa = alm_out[i].toa;
            alm_week = alm_out[i].week;
            break;
        }
    }

    if (alm_toa < 0)
    {
        alm_toa = current_bds_tow_seconds;
        alm_week = current_bds_week;
    }

    alm_toa = alignBdsAlmanacToa4096(alm_toa, &alm_week);

    for (int i = 0; i < 63; ++i)
    {
        if ((alm_out[i].valid & 1) != 0)
            continue;
        if (eph_list[i] == nullptr)
            continue;
        if ((eph_list[i]->valid & 1) == 0)
            continue;
        alm_out[i] = deriveBdsAlmanacFromEphem(eph_list[i], alm_week, alm_toa);
    }
}

int readBdsB1CEphemerisCpp(
    std::vector<ephem_t> eph_vec[MAX_SAT],
    const std::string &filename)
{
    std::ifstream infile(filename);
    if (!infile.is_open())
    {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return -1;
    }

    std::string line;
    int count = 0;

    while (std::getline(infile, line))
    {
        if (line.empty())
            continue;
        if (line[0] != 'C' && line[0] != 'B')
            continue;

        ephem_t eph{};
        datetime_t epoch{};
        double data[39]{};

        const int svid = readContentsData(line, data, 0, &epoch);
        if (svid <= 0 || svid > MAX_SAT)
            continue;

        bool record_complete = true;
        for (int row = 0; row < 7; ++row)
        {
            if (!std::getline(infile, line))
            {
                record_complete = false;
                break;
            }

            readContentsData(line, data, 3 + row * 4, nullptr);
        }

        if (!record_complete)
            break;

        eph.valid = 1;
        eph.svid = static_cast<unsigned char>(svid);

        {
            bdstime_t toc_time{};
            dateToBdt(&epoch, &toc_time);
            eph.week = toc_time.week;
            eph.toc = toc_time.milliseconds / BDS_MILLISECONDS_IN_SECOND;
        }

        eph.af0 = data[0];
        eph.af1 = data[1];
        eph.af2 = data[2];

        eph.iode = static_cast<unsigned char>(data[3]);
        eph.crs = data[4];
        eph.delta_n = data[5];
        eph.m0 = data[6];

        eph.cuc = data[7];
        eph.ecc = data[8];
        eph.cus = data[9];
        eph.sqrt_a = data[10];

        eph.toe = static_cast<int>(data[11] + 0.5);

        eph.cic = data[12];
        eph.omega0 = data[13];
        eph.cis = data[14];

        eph.i0 = data[15];
        eph.crc = data[16];
        eph.w = data[17];
        eph.omega_dot = data[18];

        eph.idot = data[19];
        eph.week = static_cast<int>(data[21]);
        eph.health = (data[24] != 0.0) ? 0x80 : 0;

        eph.tgd = data[25];

        eph.iodc = static_cast<unsigned short>(data[28]);

        eph.axis_dot = 0.0;
        eph.delta_n_dot = 0.0;

        eph.tgd_ext[0] = eph.tgd;
        eph.tgd_ext[1] = eph.tgd;
        eph.tgd_ext[2] = eph.tgd * TGD_GAMME_L5;
        eph.tgd_ext[3] = eph.tgd * TGD_GAMME_L5;
        eph.tgd_ext[4] = eph.tgd * TGD_GAMMA_E5B;

        eph.axis = eph.sqrt_a * eph.sqrt_a;
        eph.n = CGCS2000_SQRT_GM / (eph.sqrt_a * eph.axis) + eph.delta_n;
        eph.root_ecc = sqrt(1.0 - eph.ecc * eph.ecc);

        if (eph.svid <= 5 || eph.svid >= 59)
            eph.omega_delta = eph.omega_dot;
        else
            eph.omega_delta = eph.omega_dot - CGCS2000_OMEGDOTE;

        if (eph.svid <= 5 || eph.svid >= 59)
            eph.flag = 1;
        else if (eph.axis > 4e7)
            eph.flag = 2;
        else
            eph.flag = 3;

        eph_vec[eph.svid - 1].push_back(eph);
        count++;
    }

    infile.close();
    return count;
}

/**
 * \brief 为给定的观测时间找到最合适的星历数据索引
 * \param[in] obs_time 观测时间（BDT时间格式）
 * \param[in] eph 星历数据向量
 * \return 最合适的星历数据索引；如果没有找到合适的, 则返回 -1
 */
int matchEpoch(bdstime_t obs_time, vector<ephem_t> eph)
{
    int index = -1;
    double dt;
    for (unsigned i = 0; i < eph.size(); i++)
    {
        if (eph.at(i).valid == 1)
        {
            const bdstime_t toc_time = bdsTimeFromWeekSeconds(eph.at(i).week, eph.at(i).toc);
            dt = subBdsTime(obs_time, toc_time);
            if (dt >= -SECONDS_IN_HOUR && dt < SECONDS_IN_HOUR)
            {
                index = i;
                break;
            }
        }
    }
    return index;
}

/**
 * \brief 计算星历的可信时间窗口
 */
void computeBdtMinMax(
    const std::vector<ephem_t> eph_vector[MAX_SAT],
    bdstime_t *bdt_min,
    bdstime_t *bdt_max,
    datetime_t *tmin,
    datetime_t *tmax)
{
    int min_initialized = 0;
    int max_initialized = 0;
    for (int sv = 0; sv < MAX_SAT; sv++)
    {
        if (eph_vector[sv].empty())
            continue;

        const ephem_t &eph = eph_vector[sv][0];
        if (eph.valid != 1)
            continue;

        const bdstime_t toc_time = bdsTimeFromWeekSeconds(eph.week, eph.toc);
        if (!min_initialized || compareBdt(&toc_time, bdt_min) > 0)
        {
            *bdt_min = toc_time;
            min_initialized = 1;
        }
    }
    for (int sv = 0; sv < MAX_SAT; sv++)
    {
        if (eph_vector[sv].size() < 2)
            continue;

        const ephem_t &eph = eph_vector[sv][eph_vector[sv].size() - 2];
        if (eph.valid != 1)
            continue;

        const bdstime_t toc_time = bdsTimeFromWeekSeconds(eph.week, eph.toc);
        if (!max_initialized || compareBdt(&toc_time, bdt_max) > 0)
        {
            *bdt_max = toc_time;
            max_initialized = 1;
        }
    }
    if (!max_initialized && min_initialized)
    {
        *bdt_max = *bdt_min;
        max_initialized = 1;
    }
    if (min_initialized != 0)
        bdtToDate(bdt_min, tmin);
    if (max_initialized != 0)
        bdtToDate(bdt_max, tmax);
}

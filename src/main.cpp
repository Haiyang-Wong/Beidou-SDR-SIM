#ifdef _WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif
#include "bds_sdr.h"
#include <cctype>
#include <clocale>

namespace fs = std::filesystem;
static fs::path findProjectRoot(fs::path start)
{
    fs::path p = fs::absolute(start);

    while (!p.empty() && p.has_parent_path() && p != p.parent_path())
    {
        const fs::path rinex_dir = p / "rinex_files";
        const fs::path cmake_file = p / "CMakeLists.txt";

        if ((fs::exists(rinex_dir) && fs::is_directory(rinex_dir)) || fs::exists(cmake_file))
        {
            return p;
        }
        p = p.parent_path();
    }

    return fs::absolute(start);
}
static bool parseLocation(const string &s, double llh[3])
{
    return (sscanf(s.c_str(), "%lf,%lf,%lf", &llh[0], &llh[1], &llh[2]) == 3);
}

static bool parseIntArg(const char *s, int &value)
{
    char extra;
    return (sscanf(s, "%d%c", &value, &extra) == 1);
}

static bool parseDoubleArg(const char *s, double &value)
{
    char extra;
    return (sscanf(s, "%lf%c", &value, &extra) == 1);
}

static bool isValidDatetime(const datetime_t &dt)
{
    if (dt.y < 2006)
        return false;
    if (dt.m < 1 || dt.m > 12)
        return false;
    if (dt.d < 1 || dt.d > daysInMonth(dt.y, dt.m))
        return false;
    if (dt.hh < 0 || dt.hh > 23)
        return false;
    if (dt.mm < 0 || dt.mm > 59)
        return false;
    if (dt.sec < 0.0 || dt.sec >= 60.0)
        return false;

    return true;
}

static int compareDatetime(const datetime_t &a, const datetime_t &b)
{
    if (a.y != b.y)
        return a.y - b.y;
    if (a.m != b.m)
        return a.m - b.m;
    if (a.d != b.d)
        return a.d - b.d;
    if (a.hh != b.hh)
        return a.hh - b.hh;
    if (a.mm != b.mm)
        return a.mm - b.mm;
    return (a.sec > b.sec) ? 1 : (a.sec < b.sec) ? -1
                                                  : 0;
}
static int bdtUtcOffsetSeconds(const datetime_t &utc)
{
    static const datetime_t leap_effective_times[] = {
        {2009, 1, 1, 0, 0, 0.0},
        {2012, 7, 1, 0, 0, 0.0},
        {2015, 7, 1, 0, 0, 0.0},
        {2017, 1, 1, 0, 0, 0.0},
    };

    int offset = 0;
    for (const auto &leap_time : leap_effective_times)
    {
        if (compareDatetime(utc, leap_time) >= 0)
            offset++;
    }
    return offset;
}

static void normalizeBdt(bdstime_t &bdt)
{
    normalizeBdsTime(bdt);
}

static bool utcDatetimeToBdt(const datetime_t &utc, bdstime_t &bdt)
{
    if (!isValidDatetime(utc))
        return false;

    try
    {
        dateToBdt(&utc, &bdt);
    }
    catch (const std::exception &)
    {
        return false;
    }

    bdt = addBdsTime(bdt, bdtUtcOffsetSeconds(utc));
    normalizeBdt(bdt);
    return true;
}
static bool parseUtcString(const string &s, datetime_t &utc, bdstime_t &bdt)
{
    string normalized = s;
    for (char &ch : normalized)
    {
        unsigned char c = static_cast<unsigned char>(ch);
        if (ch == ',' || ch == '-' || ch == '/' || ch == ':' ||
            ch == 'T' || ch == 't' || ch == 'Z' || ch == 'z' ||
            std::isspace(c))
        {
            ch = ' ';
        }
    }

    std::istringstream iss(normalized);
    string extra;
    if (!(iss >> utc.y >> utc.m >> utc.d >> utc.hh >> utc.mm >> utc.sec))
        return false;
    if (iss >> extra)
        return false;

    return utcDatetimeToBdt(utc, bdt);
}
static bool parseUtcArgs(int argc, char *argv[], int start, int &consumed, datetime_t &utc, bdstime_t &bdt)
{
    if (start < argc && parseUtcString(argv[start], utc, bdt))
    {
        consumed = 1;
        return true;
    }

    if (start + 5 >= argc)
        return false;

    if (!parseIntArg(argv[start], utc.y) ||
        !parseIntArg(argv[start + 1], utc.m) ||
        !parseIntArg(argv[start + 2], utc.d) ||
        !parseIntArg(argv[start + 3], utc.hh) ||
        !parseIntArg(argv[start + 4], utc.mm) ||
        !parseDoubleArg(argv[start + 5], utc.sec))
    {
        return false;
    }

    if (!utcDatetimeToBdt(utc, bdt))
        return false;

    consumed = 6;
    return true;
}

static void printUsage(const char *prog)
{
    cout << "Usage: " << prog << " [options]\n\n"
         << "Options:\n"
         << "  --eph  <path>          RINEX ephemeris file path\n"
         << "                         (default: rinex_files/BRDC00IGS_R_20211700000_01D_MN.rnx)\n"
         << "  --dur  <seconds>       Simulation duration in seconds\n"
         << "                         (default: 120)\n"
         << "  --lla  <lat,lon,alt>   Static receiver position: latitude(deg),longitude(deg),altitude(m)\n"
         << "                         (default: 50,100,450)\n"
         << "  --traj <path>          Receiver trajectory CSV: time,lat,lon,alt\n"
         << "                         time is seconds from scenario start; overrides --lla\n"
         << "  --out  <path>          Output IQ file path for file generation mode\n"
         << "                         (auto-generates name if directory)\n"
         << "                         (default: auto-generated under data/output/)\n"
         << "  --rate <Hz>            Sampling rate in Hz\n"
         << "                         (default: 10.23e6)\n"
         << "  --elv  <deg>           Elevation mask in degrees\n"
         << "                         (default: 0.8)\n"
         << "  --usrp                 Stream samples to a USRP in real time; no local bin\n"
         << "                         file is written in this mode\n"
         << "  --usrp-args <args>     Optional UHD device arguments for --usrp\n"
         << "                         (default: auto-detect first UHD device)\n"
         << "  --gain <dB>            USRP TX gain for --usrp\n"
         << "                         (default: 30)\n"
         << "  --tx-freq <Hz>         USRP TX center frequency for --usrp\n"
         << "                         (default: BDS B1C 1575.42e6)\n"
         << "  --tx-prebuffer <sec>   Seconds of IQ samples to generate before USRP starts\n"
         << "                         (default: 10)\n"
         << "  --utc  <UTC time>      Start time in UTC: yyyy,mm,dd,hh,mm,ss\n"
         << "                         Also accepts: yyyy mm dd hh mm ss\n"
         << "                         (default: auto-detect from ephemeris)\n"
         << "  --help                 Show this message\n"
         << "\n"
         << "Example:\n"
         << "  " << prog << " --eph rinex_files/my.rnx --dur 37 --lla 34,108,450 --utc 2026,4,10,0,0,0\n";
}

int main(int argc, char *argv[])
{
#ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
    SetConsoleCP(CP_UTF8);
#endif
    setlocale(LC_ALL, ".UTF-8");

    option_t opt;
    sim_t sim{};
    opt.duration     = 120;
    opt.llh[0]        = 50;
    opt.llh[1]        = 100;
    opt.llh[2]        = 450;
    opt.time_overwrite = false;
    opt.g0.week       = -1;
    setBdsMillis(opt.g0, -1, 0);
    opt.samp_rate     = 10.23e6;
    opt.elv_mask      = 10;
    opt.use_trajectory = false;
    opt.traj_file[0]  = '\0';
    opt.use_usrp      = false;
    opt.tx_gain       = 30.0;
    opt.tx_freq       = CARR_FREQ;
    opt.tx_prebuffer  = 0.8;
    opt.usrp_args[0]  = '\0';

    fs::path project_root = findProjectRoot(fs::current_path());
    fs::path nav_file      = project_root / "rinex_files" / "BRDC00IGS_R_20211700000_01D_MN.rnx";
    fs::path traj_file;
    string   outfile_user;
    bool     utc_user_set = false;
    datetime_t utc_user_time{};
    for (int i = 1; i < argc; i++)
    {
        string arg = argv[i];

        if (arg == "--help" || arg == "-h")
        {
            printUsage(argv[0]);
            return 0;
        }
        else if (arg == "--eph" && i + 1 < argc)
        {
            nav_file = argv[++i];
        }
        else if (arg == "--dur" && i + 1 < argc)
        {
            opt.duration = atof(argv[++i]);
        }
        else if ((arg == "--lla" || arg == "--loc") && i + 1 < argc)
        {
            if (!parseLocation(argv[++i], opt.llh))
            {
                cerr << "Invalid --lla format. Expected: lat,lon,alt" << endl;
                return 1;
            }
        }
        else if (arg == "--traj" && i + 1 < argc)
        {
            traj_file = argv[++i];
            if (traj_file.empty())
            {
                cerr << "Invalid --traj path. Expected a CSV file path" << endl;
                return 1;
            }
            opt.use_trajectory = true;
        }
        else if (arg == "--out" && i + 1 < argc)
        {
            outfile_user = argv[++i];
        }
        else if (arg == "--rate" && i + 1 < argc)
        {
            opt.samp_rate = atof(argv[++i]);
        }
        else if (arg == "--elv" && i + 1 < argc)
        {
            opt.elv_mask = atof(argv[++i]);
        }
        else if (arg == "--usrp")
        {
            opt.use_usrp = true;
        }
        else if (arg == "--usrp-args" && i + 1 < argc)
        {
            strncpy(opt.usrp_args, argv[++i], sizeof(opt.usrp_args) - 1);
            opt.usrp_args[sizeof(opt.usrp_args) - 1] = '\0';
        }
        else if (arg == "--gain" && i + 1 < argc)
        {
            opt.tx_gain = atof(argv[++i]);
        }
        else if (arg == "--tx-freq" && i + 1 < argc)
        {
            opt.tx_freq = atof(argv[++i]);
        }
        else if (arg == "--tx-prebuffer" && i + 1 < argc)
        {
            opt.tx_prebuffer = atof(argv[++i]);
            if (opt.tx_prebuffer < 0.0)
            {
                cerr << "Invalid --tx-prebuffer. Expected a non-negative number of seconds." << endl;
                return 1;
            }
        }
        else if (arg == "--utc" && i + 1 < argc)
        {
            int consumed = 0;
            if (!parseUtcArgs(argc, argv, i + 1, consumed, utc_user_time, opt.g0))
            {
                cerr << "Invalid --utc format. Expected: yyyy,mm,dd,hh,mm,ss or yyyy mm dd hh mm ss" << endl;
                return 1;
            }
            utc_user_set = true;
            i += consumed;
        }
        else
        {
            cerr << "Unknown or incomplete option: " << arg << endl;
            printUsage(argv[0]);
            return 1;
        }
    }

    cout << "BDS SDR Simulation Started." << endl;
    time_t t = time(nullptr);
    tm tm_info{};

#ifdef _WIN32
    localtime_s(&tm_info, &t);
#else
    localtime_r(&t, &tm_info);
#endif

    char time_str[20];
    strftime(time_str, sizeof(time_str), "%m-%d_%H-%M", &tm_info);

    fs::path out_file;
    if (opt.use_usrp)
    {
        opt.out_file[0] = '\0';
        if (!outfile_user.empty())
            cout << "Note         : --out is ignored in --usrp realtime mode" << endl;
    }
    else
    {
        fs::path out_dir = project_root / "data" / "output";
        fs::create_directories(out_dir);

        if (!outfile_user.empty())
        {
            out_file = fs::path(outfile_user);
            if (fs::is_directory(out_file))
                out_file = out_file / ("bds_sdr_output_" + string(time_str) + ".bin");
        }
        else
        {
            out_file = out_dir / ("bds_sdr_output_" + string(time_str) + ".bin");
        }

        strncpy(opt.out_file, out_file.string().c_str(), sizeof(opt.out_file) - 1);
        opt.out_file[sizeof(opt.out_file) - 1] = '\0';
    }

    if (opt.use_trajectory)
    {
        fs::path candidate = traj_file;
        if (!fs::exists(candidate) && candidate.is_relative())
        {
            fs::path project_candidate = project_root / candidate;
            if (fs::exists(project_candidate))
                candidate = project_candidate;
        }

        if (!fs::exists(candidate) || !fs::is_regular_file(candidate))
        {
            cerr << "Invalid --traj path. File does not exist: " << traj_file << endl;
            return 1;
        }

        const string traj_file_string = candidate.string();
        if (traj_file_string.size() >= sizeof(opt.traj_file))
        {
            cerr << "Invalid --traj path. Path is too long." << endl;
            return 1;
        }

        strncpy(opt.traj_file, traj_file_string.c_str(), sizeof(opt.traj_file) - 1);
        opt.traj_file[sizeof(opt.traj_file) - 1] = '\0';
    }
    strncpy(opt.nav_file, nav_file.string().c_str(), sizeof(opt.nav_file) - 1);
    opt.nav_file[sizeof(opt.nav_file) - 1] = '\0';

    cout << "Project root : " << project_root << endl;
    cout << "Nav file     : " << opt.nav_file << endl;
    if (opt.use_usrp)
        cout << "Output file  : disabled in realtime mode" << endl;
    else
        cout << "Output file  : " << opt.out_file << endl;
    cout << "Duration     : " << opt.duration << " s" << endl;
    if (opt.use_trajectory)
    {
        cout << "Position     : dynamic trajectory" << endl;
        cout << "Trajectory   : " << opt.traj_file << endl;
    }
    else
        cout << "Position     : lat=" << opt.llh[0] << " lon=" << opt.llh[1] << " alt=" << opt.llh[2] << endl;
    cout << "Sample rate  : " << opt.samp_rate / 1e6 << " MHz" << endl;
    cout << "Elev mask    : " << opt.elv_mask << " deg" << endl;
    if (opt.use_usrp)
    {
        cout << "USRP TX      : enabled, freq=" << opt.tx_freq / 1e6
             << " MHz gain=" << opt.tx_gain << " dB prebuffer=" << opt.tx_prebuffer << " s";
        if (opt.usrp_args[0] == '\0')
            cout << " device=auto";
        else
            cout << " args=\"" << opt.usrp_args << "\"";
        cout << endl;
    }
    else
        cout << "USRP TX      : disabled" << endl;
    if (opt.g0.week >= 0)
    {
        if (utc_user_set)
        {
            cout << "Start time   : UTC="
                 << utc_user_time.y << "-"
                 << setw(2) << setfill('0') << utc_user_time.m << "-"
                 << setw(2) << setfill('0') << utc_user_time.d << " "
                 << setw(2) << setfill('0') << utc_user_time.hh << ":"
                 << setw(2) << setfill('0') << utc_user_time.mm << ":"
                 << fixed << setprecision(3) << setw(6) << setfill('0') << utc_user_time.sec
                 << setfill(' ') << defaultfloat << setprecision(6)
                 << " -> BDT week=" << opt.g0.week << " sec=" << ((static_cast<double>(opt.g0.milliseconds) + opt.g0.sub_milliseconds) / BDS_MILLISECONDS_IN_SECOND) << endl;
        }
        else
        {
            cout << "Start time   : week=" << opt.g0.week << " sec=" << ((static_cast<double>(opt.g0.milliseconds) + opt.g0.sub_milliseconds) / BDS_MILLISECONDS_IN_SECOND) << endl;
        }
    }
    else
        cout << "Start time   : auto-detect from ephemeris" << endl;
    sim.opt = opt;
    runBdsTask(&sim);

    return sim.status;
}

/**
 * \file main.cpp
 * \brief 主程序入口
 * \author LackWood Du
 * \date 2025-12-19
 */

#include "bds_sdr.h"
#pragma message("Compiling main.cpp")

namespace fs = std::filesystem;

// 从 start 开始向上查找项目根目录（以 rinex_files 或 CMakeLists.txt 为锚点）
static fs::path find_project_root(fs::path start)
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

    // 兜底：找不到就返回 start（至少不崩）
    return fs::absolute(start);
}

int main()
{
    cout << "BDS SDR Simulation Started." << endl;

    option_t opt;
    sim_t sim;

    // ================== 仿真参数 ==================
    opt.iduration = 37;          // 仿真时长（秒）

    // opt.llh[0] = 40.0099741;      // 纬度
    // opt.llh[1] = -105.2443329;     // 经度
    // opt.llh[2] = 1616.1;           // 高度（米）

    // opt.llh[0] = 40;      // 纬度
    // opt.llh[1] = -106;     // 经度
    // opt.llh[2] = 1000;           // 高度（米）

    opt.llh[0] = 34;      // 纬度
    opt.llh[1] = 108;     // 经度
    opt.llh[2] = 450;           // 高度（米）

    cout << "THIS IS NEW VERSION" << endl;

    opt.use_usrp = false;
    opt.timeoverwrite = false;
    opt.g0.week = -1;

    // ================== 时间处理 ==================
    time_t t = time(nullptr);
    tm tm_info{};

#ifdef _WIN32
    localtime_s(&tm_info, &t);
#else
    localtime_r(&t, &tm_info);
#endif

    char time_str[20];
    strftime(time_str, sizeof(time_str), "%m-%d_%H-%M", &tm_info);

    // ================== 路径（跨平台核心部分） ==================
    // 自动定位项目根目录（避免把 build 当成根）
    fs::path project_root = find_project_root(fs::current_path());

    // 导航文件
    // fs::path navfile = project_root / "rinex_files" / "BDS_B1C_NAV.rnx";
    // fs::path navfile = project_root / "rinex_files" / "BDS_B1C_NAV_copy.rnx";
    // fs::path navfile = project_root / "rinex_files" / "BDS_B1C_NAV_18_414M.rnx";
    fs::path navfile = project_root / "rinex_files" / "ADIS00ETH_R_20261010000_01D_CN.rnx";

    // 输出目录
    fs::path out_dir = project_root / "data" / "output";
    fs::create_directories(out_dir);

    // 输出文件
    fs::path outfile = out_dir / ("bds_sdr_output_" + string(time_str) + ".bin");

    // ================== 拷贝到 C 风格结构体 ==================
    strncpy(opt.navfile, navfile.string().c_str(), sizeof(opt.navfile) - 1);
    opt.navfile[sizeof(opt.navfile) - 1] = '\0';

    strncpy(opt.outfile, outfile.string().c_str(), sizeof(opt.outfile) - 1);
    opt.outfile[sizeof(opt.outfile) - 1] = '\0';

    // ================== 调试输出（强烈建议保留） ==================
    cout << "Project root: " << project_root << endl;
    cout << "Nav file    : " << opt.navfile << endl;
    cout << "Output file : " << opt.outfile << endl;
    cout << "Working dir : " << fs::current_path() << endl;

    // ================== 运行仿真 ==================
    sim.opt = opt;
    bds_task(&sim);

    return 0;
}
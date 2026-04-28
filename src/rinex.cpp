/**
 * \file rinex.cpp
 * \brief RINEX 文件处理相关函数实现
 * \author LackWood Du
 * \date 2025-12-19
 */

#include "bds_sdr.h"
#pragma message("Compiling rinex.cpp")

// ------------------ convertD2E ------------------
// 将 RINEX 风格的 Fortran "D" 指数转换成 "E" 科学计数法
// 例如 "-9.631020366214D-04" -> "-9.631020366214E-04"
void convertD2E(char *line)
{
    if (!line)
        return;

    for (char *p = line; *p != '\0'; ++p)
    {
        if (*p == 'D' || *p == 'd')
        {
            *p = 'E';
        }
    }
}

int readContentsData(char *str, double *data, datetime_t *time, bool read_time)
{
    int Second;
    int svid = 0;
    int length = strlen(str);

    convertD2E(str);

    if (read_time)
    {
        sscanf(str + 4, "%d %d %d %d %d %d", &(time->y), &(time->m), &(time->d), &(time->hh), &(time->mm), &Second);
        time->sec = (double)Second;
        if (str[1] == ' ')
            svid = 0;
        else
            sscanf(str + 1, "%2d", &svid);
        if (length > 24 && str[24] != ' ')
            sscanf(str + 23, "%lf", &data[0]);
        else
            data[0] = 0.0;
        if (length > 43 && str[43] != ' ')
            sscanf(str + 42, "%lf", &data[1]);
        else
            data[1] = 0.0;
        if (length > 62 && str[62] != ' ')
            sscanf(str + 61, "%lf", &data[2]);
        else
            data[2] = 0.0;
    }
    else
    {
        if (length > 5 & str[5] != ' ')
            sscanf(str + 4, "%lf", &data[0]);
        else
            data[0] = 0.0;
        if (length > 24 & str[24] != ' ')
            sscanf(str + 23, "%lf", &data[1]);
        else
            data[1] = 0.0;
        if (length > 43 && str[43] != ' ')
            sscanf(str + 42, "%lf", &data[2]);
        else
            data[2] = 0.0;
        if (length > 62 && str[62] != ' ')
            sscanf(str + 61, "%lf", &data[3]);
        else
            data[3] = 0.0;
    }

    return svid;
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
            continue; // BDS line start

        ephem_t eph{};
        datetime_t ttmp{};
        double data[4];

        char buf[200];
        memset(buf, 0, sizeof(buf));

        // -------------------- Line 1 --------------------
        strncpy(buf, line.c_str(), sizeof(buf) - 1);
        int svid = readContentsData(buf, data, &ttmp, true);

        eph.PRN = svid;
        eph.svid = svid;

        eph.af0 = data[0];
        eph.af1 = data[1];
        eph.af2 = data[2];

        // convert UTC to BDT toc
        date2bdt(&ttmp, &eph.toc);

        // -------------------- Line 2 --------------------
        if (!std::getline(infile, line))
            break;
        strncpy(buf, line.c_str(), sizeof(buf) - 1);
        readContentsData(buf, data, nullptr, false);

        eph.iode = (int)data[0];
        eph.crs = data[1];
        eph.deltan = data[2];
        eph.m0 = data[3];

        // -------------------- Line 3 --------------------
        if (!std::getline(infile, line))
            break;
        strncpy(buf, line.c_str(), sizeof(buf) - 1);
        readContentsData(buf, data, nullptr, false);

        eph.cuc = data[0];
        eph.ecc = data[1];
        eph.cus = data[2];
        eph.sqrta = data[3];

        // -------------------- Line 4 --------------------
        if (!std::getline(infile, line))
            break;
        strncpy(buf, line.c_str(), sizeof(buf) - 1);
        readContentsData(buf, data, nullptr, false);

        eph.toe.sec = (int)(data[0] + 0.5);
        eph.cic = data[1];
        eph.omg0 = data[2];
        eph.cis = data[3];

        // toe.week 无法从 B1C 获取，先置 0
        eph.toe.week = 0;

        // -------------------- Line 5 --------------------
        if (!std::getline(infile, line))
            break;
        strncpy(buf, line.c_str(), sizeof(buf) - 1);
        readContentsData(buf, data, nullptr, false);

        eph.inc0 = data[0];
        eph.crc = data[1];
        eph.aop = data[2];
        eph.omgdot = data[3];

        // -------------------- Line 6 --------------------
        if (!std::getline(infile, line))
            break;
        strncpy(buf, line.c_str(), sizeof(buf) - 1);
        readContentsData(buf, data, nullptr, false);

        eph.idot = data[0];
        // data[1] = spare
        eph.toe.week = (int)data[2]; // toe week
        // data[3] = spare

        // -------------------- Line 7 --------------------
        if (!std::getline(infile, line))
            break;
        strncpy(buf, line.c_str(), sizeof(buf) - 1);
        readContentsData(buf, data, nullptr, false);

        eph.ura = data[0];
        eph.svhlth = data[1];
        eph.tgd1 = data[2];
        eph.tgd2 = data[3];

        // -------------------- Line 8 --------------------
        if (!std::getline(infile, line))
            break;
        strncpy(buf, line.c_str(), sizeof(buf) - 1);
        readContentsData(buf, data, nullptr, false);

        eph.gps_time = data[0];
        eph.IODC = (int)data[1];
        // data[2], data[3] = spare

        // ---------------- Derived parameters ----------------
        eph.A = eph.sqrta * eph.sqrta;
        eph.n = WGS_SQRT_GM / (eph.A * eph.sqrta) + eph.deltan;
        eph.sq1e2 = sqrt(1.0 - eph.ecc * eph.ecc);
        eph.omg_t = eph.omg0 - OMEGA_EARTH * eph.toe.sec;
        eph.omgkdot = eph.omgdot - OMEGA_EARTH;
        eph.vflg = 1;

        eph_vec[eph.PRN - 1].push_back(eph);
        count++;
    }

    infile.close();
    return count;
}

// ------------------ printEphVec ------------------
// 打印卫星星历信息
void printEphVec(const std::vector<ephem_t> eph_vec[MAX_SAT])
{
    for (int prn = 0; prn < MAX_SAT; prn++)
    {
        if (eph_vec[prn].empty())
            continue;

        std::cout << "========== PRN " << prn + 1 << " ==========\n";

        for (size_t i = 0; i < eph_vec[prn].size(); i++)
        {
            const ephem_t &e = eph_vec[prn][i];

            std::cout << "\n----- Ephemeris #" << i + 1 << " -----\n";

            // -------- Line 1 --------
            std::cout << "Line1 PRN: " << e.PRN << "\n";
            std::cout << "Line1 af0: " << e.af0 << "\n";
            std::cout << "Line1 af1: " << e.af1 << "\n";
            std::cout << "Line1 af2: " << e.af2 << "\n";
            std::cout << "Line1 toc (week,sec): "
                      << e.toc.week << ", " << e.toc.sec << "\n";

            // -------- Line 2 --------
            std::cout << "Line2 iode: " << e.iode << "\n";
            std::cout << "Line2 crs: " << e.crs << "\n";
            std::cout << "Line2 deltan: " << e.deltan << "\n";
            std::cout << "Line2 m0: " << e.m0 << "\n";

            // -------- Line 3 --------
            std::cout << "Line3 cuc: " << e.cuc << "\n";
            std::cout << "Line3 ecc: " << e.ecc << "\n";
            std::cout << "Line3 cus: " << e.cus << "\n";
            std::cout << "Line3 sqrtA: " << e.sqrta << "\n";

            // -------- Line 4 --------
            std::cout << "Line4 toe.sec: " << e.toe.sec << "\n";
            std::cout << "Line4 cic: " << e.cic << "\n";
            std::cout << "Line4 omg0: " << e.omg0 << "\n";
            std::cout << "Line4 cis: " << e.cis << "\n";

            // -------- Line 5 --------
            std::cout << "Line5 inc0: " << e.inc0 << "\n";
            std::cout << "Line5 crc: " << e.crc << "\n";
            std::cout << "Line5 aop: " << e.aop << "\n";
            std::cout << "Line5 omgdot: " << e.omgdot << "\n";

            // -------- Line 6 --------
            std::cout << "Line6 idot: " << e.idot << "\n";

            // -------- Line 7 --------
            std::cout << "Line7 ura: " << e.ura << "\n";
            std::cout << "Line7 svhealth: " << e.svhlth << "\n";
            std::cout << "Line7 tgd1: " << e.tgd1 << "\n";
            std::cout << "Line7 tgd2: " << e.tgd2 << "\n";

            // -------- Line 8 --------
            std::cout << "Line8 gps_time: " << e.gps_time << "\n";
            std::cout << "Line8 IODC: " << e.IODC << "\n";

            // -------- Derived parameters --------
            std::cout << "Derived A: " << e.A << "\n";
            std::cout << "Derived n: " << e.n << "\n";
            std::cout << "Derived sq1e2: " << e.sq1e2 << "\n";
            std::cout << "Derived omg_t: " << e.omg_t << "\n";
            std::cout << "Derived omgkdot: " << e.omgkdot << "\n";
            std::cout << "Derived vflg: " << e.vflg << "\n";
        }
    }
}

/**
 * \brief 为给定的观测时间找到最合适的星历数据索引
 * \param[in] obsTime 观测时间（BDT时间格式）
 * \param[in] eph 星历数据向量
 * \return 最合适的星历数据索引；如果没有找到合适的, 则返回 -1
 */
int epoch_matcher(bdstime_t obsTime, vector<ephem_t> eph)
{
    int index = -1;
    double dt;
    // 遍历所有的星历数据
    for (unsigned i = 0; i < eph.size(); i++)
    {
        // 判断当前星历是否有效
        if (eph.at(i).vflg == 1)
        {
            // 计算观测时间和星历参考时间之间的差值, BDS星历更新周期为一个小时
            dt = subBdsTime(obsTime, eph.at(i).toc);

            // 如果时间差在±1小时范围内，认为找到了合适的星历
            if (dt >= -SECONDS_IN_HOUR && dt < SECONDS_IN_HOUR)
            {
                // cout << "EM: " << dt << " - " << eph.at(i).svid + 1 << endl;
                index = i;
                break;
            }
        }
    }
    // Index of most appropriate Nav vector
    // 返回合适的星历数据索引
    return index;
}

/**
 * \brief 计算星历的可信时间窗口
 */
void compute_bdt_min_max(
    const std::vector<ephem_t> eph_vector[MAX_SAT],
    bdstime_t *bdt_min,
    bdstime_t *bdt_max,
    datetime_t *tmin,
    datetime_t *tmax)
{
    int sv;
    int min_initialized = 0;
    int max_initialized = 0;
    ephem_t eph;

    /* ---------- 计算 bdt_min ---------- */
    for (sv = 0; sv < MAX_SAT; sv++)
    {
        if (eph_vector[sv].empty())
            continue;

        eph = eph_vector[sv][0]; // 第一条星历（假定已按时间排序）

        if (eph.vflg == 1)
        {
            if (!min_initialized)
            {
                *bdt_min = eph.toc;
                min_initialized = 1;
            }
            else if (bdt_cmp(&eph.toc, bdt_min) > 0)
            {
                *bdt_min = eph.toc;
            }
        }
    }

    /* ---------- 计算 bdt_max ---------- */
    for (sv = 0; sv < MAX_SAT; sv++)
    {
        if (eph_vector[sv].size() < 2)
            continue;

        eph = eph_vector[sv][eph_vector[sv].size() - 2]; // 倒数第二条

        if (eph.vflg == 1)
        {
            if (!max_initialized)
            {
                *bdt_max = eph.toc;
                max_initialized = 1;
            }
            else if (bdt_cmp(&eph.toc, bdt_max) > 0)
            {
                *bdt_max = eph.toc;
            }
        }
    }
    // 转换为日期时间结构体
    if (min_initialized != 0)
        bdt2date(bdt_min, tmin);
    if (max_initialized != 0)
        bdt2date(bdt_max, tmax);
}
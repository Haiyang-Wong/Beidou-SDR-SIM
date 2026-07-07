#include "bds_sdr.h"

double subBdsTime(bdstime_t g1, bdstime_t g0)
{
    const long long week_ms = static_cast<long long>(g1.week - g0.week) * BDS_MILLISECONDS_IN_WEEK;
    const long long millisecond_delta = week_ms +
                                        static_cast<long long>(g1.milliseconds) -
                                        static_cast<long long>(g0.milliseconds);
    const double sub_millisecond_delta = g1.sub_milliseconds - g0.sub_milliseconds;

    return (static_cast<double>(millisecond_delta) + sub_millisecond_delta) /
           BDS_MILLISECONDS_IN_SECOND;
}
bool isLeapYear(int year)
{
    return (year % 4 == 0 && year % 100 != 0) || (year % 400 == 0);
}
int daysInMonth(int year, int month)
{
    if (month != 2)
        return mdays[month - 1];
    return mdays[1] + (isLeapYear(year) ? 1 : 0);
}
void dateToBdt(const datetime_t *dt, bdstime_t *bdt)
{
    int ref_year = 2006;
    int ref_month = 1;
    int ref_day = 1;
    int days = 0;
    for (int y = ref_year; y < dt->y; ++y)
    {
        days += isLeapYear(y) ? 366 : 365;
    }
    for (int m = 1; m < dt->m; ++m)
    {
        days += daysInMonth(dt->y, m);
    }
    days += dt->d - ref_day;

    if (days < 0)
    {
        throw std::runtime_error("datetime earlier than BDS epoch (2006-01-01)");
    }

    const int week = days / 7;
    const int day_of_week = days % 7;
    const double second_ms = dt->sec * BDS_MILLISECONDS_IN_SECOND;
    const long long whole_second_ms = static_cast<long long>(std::floor(second_ms));
    const double sub_ms = second_ms - static_cast<double>(whole_second_ms);
    const long long milliseconds_of_week =
        static_cast<long long>(day_of_week) * BDS_MILLISECONDS_IN_DAY +
        static_cast<long long>(dt->hh) * BDS_MILLISECONDS_IN_HOUR +
        static_cast<long long>(dt->mm) * BDS_MILLISECONDS_IN_MINUTE +
        whole_second_ms;

    setBdsMillis(*bdt, week, milliseconds_of_week, sub_ms);
}


void bdtToDate(const bdstime_t *bdt, datetime_t *date)
{
    int days, sec_day;
    int year = 2006, month = 1;
    int d;

    bdstime_t normalized = *bdt;
    normalizeBdsTime(normalized);

    /* BDT 起点：2006-01-01 */
    days = normalized.week * 7 + normalized.milliseconds / BDS_MILLISECONDS_IN_DAY;
    const int ms_day = normalized.milliseconds % BDS_MILLISECONDS_IN_DAY;
    sec_day = ms_day / BDS_MILLISECONDS_IN_SECOND;

    /* 计算年 */
    while (1)
    {
        int ydays = isLeapYear(year) ? 366 : 365;
        if (days >= ydays)
        {
            days -= ydays;
            year++;
        }
        else
        {
            break;
        }
    }

    /* 计算月 */
    for (month = 0; month < 12; month++)
    {
        d = mdays[month];
        if (month == 1 && isLeapYear(year))
            d++;

        if (days >= d)
        {
            days -= d;
        }
        else
        {
            break;
        }
    }

    /* 填充 date 结构 */
    date->y = year;
    date->m = month + 1;
    date->d = days + 1;

    date->hh = sec_day / 3600;
    sec_day %= 3600;
    date->mm = sec_day / 60;
    const int second_in_minute = sec_day % 60;
    const int millisecond_in_second = ms_day % BDS_MILLISECONDS_IN_SECOND;
    date->sec = static_cast<double>(second_in_minute) +
                (static_cast<double>(millisecond_in_second) + normalized.sub_milliseconds) /
                    BDS_MILLISECONDS_IN_SECOND;
}


int compareBdt(const bdstime_t *t1, const bdstime_t *t2)
{
    bdstime_t lhs = *t1;
    bdstime_t rhs = *t2;
    normalizeBdsTime(lhs);
    normalizeBdsTime(rhs);

    if (lhs.week != rhs.week)
        return lhs.week - rhs.week;
    if (lhs.milliseconds != rhs.milliseconds)
        return (lhs.milliseconds > rhs.milliseconds) ? 1 : -1;
    return (lhs.sub_milliseconds > rhs.sub_milliseconds) ? 1 :
           (lhs.sub_milliseconds < rhs.sub_milliseconds) ? -1 :
                                                          0;
}


/**
 * \brief 为一个给定的BDS时间点增加一个指定的时间增量，得到一个新的BDS时间点
 * \param[in] g 是输入的 GAL 日期时间结构体
 * \param[in] dt 是要增加的时间增量，单位为秒，可以是正数或负数
 * \return 返回增加时间后的新的 GAL 日期时间结构体
 */
bdstime_t addBdsTime(bdstime_t g, double dt)
{
    const double delta_ms = dt * BDS_MILLISECONDS_IN_SECOND;
    const double whole_delta_ms = std::floor(delta_ms);
    const long long milliseconds = static_cast<long long>(g.milliseconds) +
                                   static_cast<long long>(whole_delta_ms);
    const double sub_milliseconds = g.sub_milliseconds + delta_ms - whole_delta_ms;

    setBdsMillis(g, g.week, milliseconds, sub_milliseconds);
    return g;
}


/**
 * \brief 解决的核心问题是：仿真从什么时候开始（信号从什么时候开始生成）？这个时间是否在星历有效时间内？是否需要把星历整体平移到该时间
 * \param[in, out] b0 用户指定的仿真场景起始时间，如果用户未指定则使用星历的最早时间
 * \param[in] bmin 星历的最早时间
 * \param[in] bmax 星历的最晚时间
 * \param[out] t0 对应 b0 的日期时间结构体（公历时间）
 * \param[out] tmin 对应 bmin 的日期时间结构体（公历时间）
 * \param[out] tmax 对应 bmax 的日期时间结构体（公历时间）
 * \param[in] time_overwrite 如果为 true，表示用户希望强制对齐星历时间到指定的 b0 时间
 * \param[in] neph 星历数据的数量
 * \param[in, out] eph1 星历数据数组
 */
void setScenarioStartTimeBds(bdstime_t *b0,
                                 bdstime_t bmin,
                                 bdstime_t bmax,
                                 datetime_t *t0,
                                 datetime_t *tmin,
                                 datetime_t *tmax,
                                 bool time_overwrite,
                                 int neph,
                                 vector<ephem_t> eph1[MAX_SAT])
{
    if (b0->week >= 0) // 用户已设置场景起始时间
    {
        if (time_overwrite == true) // 重新对齐星历时间
        {
            bdstime_t btmp;
            double dsec;

            /* 对齐到 2 小时边界（与星历参考时间一致） */
            btmp.week = b0->week;
            btmp.milliseconds = (b0->milliseconds / BDS_MILLISECONDS_IN_HOUR) * BDS_MILLISECONDS_IN_HOUR;
            btmp.sub_milliseconds = 0.0;

            /* 计算相对于星历最早时间的平移量 */
            dsec = subBdsTime(btmp, bmin);

            /* 对所有卫星、所有星历执行时间平移 */
            for (int sv = 0; sv < MAX_SAT; sv++)
            {
                for (int i = 0; i < neph; i++)
                {
                    if (eph1[i][sv].valid == 1)
                    {
                        /* TOC 对齐 */
                        btmp = bdsTimeFromWeekSeconds(eph1[i][sv].week, eph1[i][sv].toc);
                        btmp = addBdsTime(btmp, dsec);
                        eph1[i][sv].week = btmp.week;
                        eph1[i][sv].toc = btmp.milliseconds / BDS_MILLISECONDS_IN_SECOND;

                        /* TOE 对齐 */
                        btmp = bdsTimeFromWeekSeconds(eph1[i][sv].week, eph1[i][sv].toe);
                        btmp = addBdsTime(btmp, dsec);
                        eph1[i][sv].week = btmp.week;
                        eph1[i][sv].toe = btmp.milliseconds / BDS_MILLISECONDS_IN_SECOND;
                    }
                }
            }
        }
        else // 仅检查时间是否合法
        {
            if (subBdsTime(*b0, bmin) < 0.0 ||
                subBdsTime(bmax, *b0) < 0.0)
            {
                datetime_t tl;
                bdtToDate(b0, &tl);

                fprintf(stderr,
                        "Start = %4d/%02d/%02d,%02d:%02d:%02.0f (%d:%.0f)\n",
                        tl.y, tl.m, tl.d, tl.hh, tl.mm, tl.sec,
                        b0->week, (static_cast<double>(b0->milliseconds) + b0->sub_milliseconds) / BDS_MILLISECONDS_IN_SECOND);

                fprintf(stderr, "ERROR: Invalid start time.\n");

                fprintf(stderr,
                        "tmin = %4d/%02d/%02d,%02d:%02d:%02.0f (%d:%.0f)\n",
                        tmin->y, tmin->m, tmin->d,
                        tmin->hh, tmin->mm, tmin->sec,
                        bmin.week, (static_cast<double>(bmin.milliseconds) + bmin.sub_milliseconds) / BDS_MILLISECONDS_IN_SECOND);

                fprintf(stderr,
                        "tmax = %4d/%02d/%02d,%02d:%02d:%02.0f (%d:%.0f)\n",
                        tmax->y, tmax->m, tmax->d,
                        tmax->hh, tmax->mm, tmax->sec,
                        bmax.week, (static_cast<double>(bmax.milliseconds) + bmax.sub_milliseconds) / BDS_MILLISECONDS_IN_SECOND);

                exit(1);
            }
        }
    }
    else // 用户未设置起始时间，使用最早星历时间
    {
        b0->week = bmin.week;
        *b0 = bmin;

        t0 = tmin;
    }

    datetime_t tl;

    bdtToDate(b0, &tl);

    cout << "------------------ Ephemeris valid time range -------------------" << endl;
    cout << "tmin = "
         << setw(4) << setfill('0') << tmin->y << "/"
         << setw(2) << setfill('0') << tmin->m << "/"
         << setw(2) << setfill('0') << tmin->d << ","
         << setw(2) << setfill('0') << tmin->hh << ":"
         << setw(2) << setfill('0') << tmin->mm << ":"
         << setw(2) << setfill('0') << fixed << setprecision(0) << tmin->sec
         << " (" << bmin.week << ":" << defaultfloat << setprecision(6) << ((static_cast<double>(bmin.milliseconds) + bmin.sub_milliseconds) / BDS_MILLISECONDS_IN_SECOND) << ")" << endl;
    cout << "tmax = "
         << setw(4) << setfill('0') << tmax->y << "/"
         << setw(2) << setfill('0') << tmax->m << "/"
         << setw(2) << setfill('0') << tmax->d << ","
         << setw(2) << setfill('0') << tmax->hh << ":"
         << setw(2) << setfill('0') << tmax->mm << ":"
         << setw(2) << setfill('0') << fixed << setprecision(0) << tmax->sec
         << " (" << bmax.week << ":" << defaultfloat << setprecision(6) << ((static_cast<double>(bmax.milliseconds) + bmax.sub_milliseconds) / BDS_MILLISECONDS_IN_SECOND) << ")" << endl;
    cout << "------------------ Scenario start time -------------------" << endl;
    cout << "Start = "
         << setw(4) << setfill('0') << tl.y << "/"
         << setw(2) << setfill('0') << tl.m << "/"
         << setw(2) << setfill('0') << tl.d << ","
         << setw(2) << setfill('0') << tl.hh << ":"
         << setw(2) << setfill('0') << tl.mm << ":"
         << setw(2) << setfill('0') << fixed << setprecision(0) << tl.sec
         << " (" << b0->week << ":" << defaultfloat << setprecision(6) << ((static_cast<double>(b0->milliseconds) + b0->sub_milliseconds) / BDS_MILLISECONDS_IN_SECOND) << ")"
         << setfill(' ') << endl;
}

/**
 * \brief 根据接收时间和伪距，计算卫星信号发射时间
 * \param[in] grx 接收机BDT时间
 * \param[in] pseudorange 伪距（米）
 * \return 卫星发射时的BDT时间
 */
bdstime_t computeSatelliteTxTime(bdstime_t grx, double pseudorange)
{
    const double delay = pseudorange / SPEED_OF_LIGHT;
    return addBdsTime(grx, -delay);
}

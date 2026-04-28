/**
 * @file bds_sdr.cpp
 * @brief BDS SDR 时间处理相关函数实现
 * @author LackWood Du
 * @date 2025-12-23
 */

#include "bds_sdr.h"
#pragma message("Compiling gnss_time.cpp")

double subBdsTime(bdstime_t g1, bdstime_t g0)
{
    double dt;

    dt = g1.sec - g0.sec;
    dt += (double)(g1.week - g0.week) * SECONDS_IN_WEEK;

    return (dt);
}

// ------------------ date2bdt ------------------
// 判断闰年
bool isLeapYear(int year)
{
    return (year % 4 == 0 && year % 100 != 0) || (year % 400 == 0);
}

// 返回当年当月天数
int daysInMonth(int year, int month)
{
    if (month != 2)
        return mdays[month - 1];
    return mdays[1] + (isLeapYear(year) ? 1 : 0);
}

// 精确版 datetime_t -> bdstime_t
void date2bdt(const datetime_t *dt, bdstime_t *bdt)
{
    // 北斗时间参考：BDS周起点 2006-01-01 00:00:00
    int ref_year = 2006;
    int ref_month = 1;
    int ref_day = 1;

    // 计算从参考日到目标日期的天数
    int days = 0;

    // 年
    for (int y = ref_year; y < dt->y; ++y)
    {
        days += isLeapYear(y) ? 366 : 365;
    }

    // 月
    for (int m = 1; m < dt->m; ++m)
    {
        days += daysInMonth(dt->y, m);
    }

    // 日
    days += dt->d - ref_day;

    if (days < 0)
    {
        throw std::runtime_error("datetime earlier than BDS epoch (2006-01-01)");
    }

    // 计算周号与周内秒
    bdt->week = days / 7;
    int day_of_week = days % 7;

    bdt->sec = day_of_week * 86400.0 + dt->hh * 3600.0 + dt->mm * 60.0 + dt->sec;
}

void bdt2date(const bdstime_t *bdt, datetime_t *date)
{
    int days, sec_day;
    int year = 2006, month = 1;
    int d;

    /* BDT 起点：2006-01-01 */
    days = bdt->week * 7 + (int)(bdt->sec / 86400.0);
    sec_day = (int)(bdt->sec - days % 7 * 86400);

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
    date->sec = sec_day % 60 + (bdt->sec - floor(bdt->sec));
}

int bdt_cmp(const bdstime_t *t1, const bdstime_t *t2)
{
    if (t1->week != t2->week)
        return t1->week - t2->week;
    return (t1->sec > t2->sec) ? 1 : (t1->sec < t2->sec) ? -1
                                                         : 0;
}

/**
 * \brief 为一个给定的GPS时间点增加一个指定的时间增量，得到一个新的GPS时间点
 * \param[in] g 是输入的 GAL 日期时间结构体
 * \param[in] dt 是要增加的时间增量，单位为秒，可以是正数或负数
 * \return 返回增加时间后的新的 GAL 日期时间结构体
 */
bdstime_t incBdsTime(bdstime_t g, double dt)
{
    g.sec = g.sec + dt;
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
 * \param[in] timeoverwrite 如果为 true，表示用户希望强制对齐星历时间到指定的 b0 时间
 * \param[in] neph 星历数据的数量
 * \param[in, out] eph1 星历数据数组
 */
void set_scenario_start_time_bds(bdstime_t *b0,
                                 bdstime_t bmin,
                                 bdstime_t bmax,
                                 datetime_t *t0,
                                 datetime_t *tmin,
                                 datetime_t *tmax,
                                 bool timeoverwrite,
                                 int neph,
                                 vector<ephem_t> eph1[MAX_SAT])
{
    if (b0->week >= 0) // 用户已设置场景起始时间
    {
        // 用户给了起始时间，但是这个时间不一定和星历参考时间一致，因此需要整体平移星历，使其看起来像是从该时间开始生成的
        if (timeoverwrite == true) // 重新对齐星历时间
        {
            bdstime_t btmp;
            datetime_t ttmp;
            double dsec;

            /* 对齐到 2 小时边界（与星历参考时间一致） */
            btmp.week = b0->week;
            //  秒 = 1 小时，以1小时为一个周期，当前时刻所在的周期起点
            btmp.sec = (double)(((int)(b0->sec)) / 3600) * 3600.0;

            /* 计算相对于星历最早时间的平移量 */
            // dsec = btmp - bmin，单位是秒
            dsec = subBdsTime(btmp, bmin);

            /* 对所有卫星、所有星历执行时间平移 */
            for (int sv = 0; sv < MAX_SAT; sv++)
            {
                for (int i = 0; i < neph; i++)
                {
                    if (eph1[i][sv].vflg == 1)
                    {
                        /* TOC 对齐 */
                        btmp = incBdsTime(eph1[i][sv].toc, dsec);
                        bdt2date(&btmp, &ttmp);
                        eph1[i][sv].toc = btmp;
                        eph1[i][sv].t = ttmp;

                        /* TOE 对齐 */
                        btmp = incBdsTime(eph1[i][sv].toe, dsec);
                        eph1[i][sv].toe = btmp;
                    }
                }
            }
        }
        else // 仅检查时间是否合法
        {
            // 如果不合法则报错退出
            // 这里严格将仿真起始时间限定在bmin和bmax之间，是为了保证仿真过程中所有时刻都有有效星历
            // 而非仅仅考虑某一条星历的有效时间是在(toe - 1, toe + 1)小时范围内
            if (subBdsTime(*b0, bmin) < 0.0 ||
                subBdsTime(bmax, *b0) < 0.0)
            {
                datetime_t tl;
                bdt2date(b0, &tl);

                fprintf(stderr,
                        "Start = %4d/%02d/%02d,%02d:%02d:%02.0f (%d:%.0f)\n",
                        tl.y, tl.m, tl.d, tl.hh, tl.mm, tl.sec,
                        b0->week, b0->sec);

                fprintf(stderr, "ERROR: Invalid start time.\n");

                fprintf(stderr,
                        "tmin = %4d/%02d/%02d,%02d:%02d:%02.0f (%d:%.0f)\n",
                        tmin->y, tmin->m, tmin->d,
                        tmin->hh, tmin->mm, tmin->sec,
                        bmin.week, bmin.sec);

                fprintf(stderr,
                        "tmax = %4d/%02d/%02d,%02d:%02d:%02.0f (%d:%.0f)\n",
                        tmax->y, tmax->m, tmax->d,
                        tmax->hh, tmax->mm, tmax->sec,
                        bmax.week, bmax.sec);

                exit(1);
            }
        }
    }
    else // 用户未设置起始时间，使用最早星历时间
    {
        b0->week = bmin.week;
        b0->sec = bmin.sec;

        t0 = tmin;
    }

    datetime_t tl;

    bdt2date(b0, &tl);

    cout << "------------------星历数据的有效时间范围：-------------------" << endl;
    fprintf(stderr, "\ntmin = %4d/%02d/%02d,%02d:%02d:%02.0f (%d:%.0f)\n", tmin->y,
            tmin->m, tmin->d, tmin->hh, tmin->mm, tmin->sec, bmin.week, bmin.sec);
    fprintf(stderr, "tmax = %4d/%02d/%02d,%02d:%02d:%02.0f (%d:%.0f)\n", tmax->y,
            tmax->m, tmax->d, tmax->hh, tmax->mm, tmax->sec, bmax.week, bmax.sec);
    cout << "------------------仿真场景起始时间设置为：-------------------" << endl;
    fprintf(stderr, "Start = %4d/%02d/%02d,%02d:%02d:%02.0f (%d:%.0f)\n ",
            tl.y, tl.m, tl.d, tl.hh, tl.mm, tl.sec, b0->week, b0->sec);
}

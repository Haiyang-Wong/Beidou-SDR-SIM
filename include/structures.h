#pragma once

#include "bds_sdr.h"

struct datetime_t
{
    int y;
    int m;
    int d;
    int hh;
    int mm;
    double sec;
};

struct bdstime_t
{
    int week;
    int milliseconds;
    double sub_milliseconds;
};

constexpr int BDS_MILLISECONDS_IN_SECOND = 1000;
constexpr int BDS_MILLISECONDS_IN_MINUTE = 60 * BDS_MILLISECONDS_IN_SECOND;
constexpr int BDS_MILLISECONDS_IN_HOUR = 60 * BDS_MILLISECONDS_IN_MINUTE;
constexpr int BDS_MILLISECONDS_IN_DAY = 24 * BDS_MILLISECONDS_IN_HOUR;
constexpr int BDS_MILLISECONDS_IN_WEEK = 7 * BDS_MILLISECONDS_IN_DAY;
constexpr int BDS_MILLISECONDS_IN_HALF_WEEK = BDS_MILLISECONDS_IN_WEEK / 2;

inline void normalizeBdsTime(bdstime_t &time)
{
    long long milliseconds = time.milliseconds;
    double sub_milliseconds = time.sub_milliseconds;

    if (sub_milliseconds >= 1.0 || sub_milliseconds < 0.0)
    {
        const double whole_sub_ms = std::floor(sub_milliseconds);
        milliseconds += static_cast<long long>(whole_sub_ms);
        sub_milliseconds -= whole_sub_ms;
    }

    if (milliseconds >= BDS_MILLISECONDS_IN_WEEK || milliseconds < 0)
    {
        long long week_delta = milliseconds / BDS_MILLISECONDS_IN_WEEK;
        if (milliseconds < 0 && (milliseconds % BDS_MILLISECONDS_IN_WEEK) != 0)
            week_delta--;

        time.week += static_cast<int>(week_delta);
        milliseconds -= week_delta * BDS_MILLISECONDS_IN_WEEK;
    }

    time.milliseconds = static_cast<int>(milliseconds);
    time.sub_milliseconds = sub_milliseconds;
}

inline void setBdsMillis(bdstime_t &time, int week, long long milliseconds_of_week, double sub_milliseconds = 0.0)
{
    time.week = week;
    time.milliseconds = 0;
    time.sub_milliseconds = 0.0;

    long long milliseconds = milliseconds_of_week;
    if (sub_milliseconds >= 1.0 || sub_milliseconds < 0.0)
    {
        const double whole_sub_ms = std::floor(sub_milliseconds);
        milliseconds += static_cast<long long>(whole_sub_ms);
        sub_milliseconds -= whole_sub_ms;
    }

    if (milliseconds >= BDS_MILLISECONDS_IN_WEEK || milliseconds < 0)
    {
        long long week_delta = milliseconds / BDS_MILLISECONDS_IN_WEEK;
        if (milliseconds < 0 && (milliseconds % BDS_MILLISECONDS_IN_WEEK) != 0)
            week_delta--;

        time.week += static_cast<int>(week_delta);
        milliseconds -= week_delta * BDS_MILLISECONDS_IN_WEEK;
    }

    time.milliseconds = static_cast<int>(milliseconds);
    time.sub_milliseconds = sub_milliseconds;
}

inline bdstime_t bdsTimeFromWeekSeconds(int week, int seconds_of_week)
{
    bdstime_t time{};
    setBdsMillis(time, week, static_cast<long long>(seconds_of_week) * BDS_MILLISECONDS_IN_SECOND);
    return time;
}

struct ephem_t
{
    unsigned short iodc;
    unsigned char iode;
    unsigned char svid;
    unsigned char valid;
    unsigned short flag;
    unsigned short health;

    int toe;
    int toc;
    int toa;
    int week;

    double m0;
    double delta_n;
    double delta_n_dot;
    double ecc;
    double sqrt_a;
    double axis_dot;
    double omega0;
    double i0;
    double w;
    double omega_dot;
    double idot;
    double cuc;
    double cus;
    double crc;
    double crs;
    double cic;
    double cis;

    double af0;
    double af1;
    double af2;
    double tgd;
    double tgd_ext[5];

    double axis;
    double n;
    double root_ecc;
    double omega_delta;
};

struct range_t
{
    bdstime_t g;
    double range;
    double d;
    double azel[2];
};

struct channel_t
{
    int prn;
    short *ca_b1c_data;
    short *ca_b1c_pilot11;
    short *ca_b1c_pilot61;
    double *weighted_b1c_data;
    double *weighted_b1c_pilot11;
    double *weighted_b1c_pilot61;
    short *b1c_pilot_sub_code;
    double f_carr;
    double f_code;
    short *nav_bit;
    double carr_phase;
    int ibit;
    double azel[2];
    range_t rho0;
    bool set_code_phase;
    double code_phase;
};

struct option_t
{
    char nav_file[MAX_CHAR];
    char out_file[MAX_CHAR];
    double duration;
    bdstime_t g0;
    double llh[3];
    int time_overwrite;
    double samp_rate;
    double elv_mask;
    bool use_trajectory;
    char traj_file[MAX_CHAR];
    bool use_usrp;
    double tx_gain;
    double tx_freq;
    double tx_prebuffer;
    char usrp_args[MAX_CHAR];
};

struct sim_t
{
    option_t opt;
    int status;
};

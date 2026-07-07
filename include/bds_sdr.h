#pragma once

#include <iostream>
using namespace std;
#include <vector>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <cstring>
#include <stdexcept>
#include <algorithm>
#include <cstdint>
#include <time.h>
#include <filesystem>
#include "constants.h"
#include "structures.h"

struct b1c_nav_time_fields_t
{
    double tx_ms_raw;
    double tx_ms_normalized;
    double hour_start_ms;
    double frame_start_ms;

    double tx_sow_raw;
    double tx_sow_normalized;
    double hour_start_sow;
    double frame_ratio;
    int floor_result;
    int how;
    int soh;
    int tow_from_how_soh;
    double frame_start_sow;
};

// geodesy函数声明
bool isGEO(const int &prn);
void satPos(ephem_t eph, bdstime_t g, double *pos, double *vel, double *clk);
void subVector(double *y, const double *x1, const double *x2);
double normVector(const double *x);
void xyz2llh(const double *xyz, double *llh);
void llh2xyz(const double *llh, double *xyz);
void ltcmat(const double *llh, double t[3][3]);
void ecef2neu(const double *xyz, double t[3][3], double *neu);
void neu2azel(double *azel, const double *neu);
int checkSatVisibility(ephem_t eph, bdstime_t g, double *xyz, double *azel, int prn, double elv_mask);

// bds_sig函数声明
void computeRange(range_t *rho, ephem_t eph, bdstime_t g, double xyz[], int prn);
void computeCodePhase(channel_t *chan, range_t rho1, double dt, bdstime_t grx);
void generateB1CDataCode(short *ca, int prn);
void generateB1CPilotCode(short *ca11, short *ca61, int prn);
void hexToB1CCa(std::vector<short> &tmp_ca, int prn, int flag);
void generateB1CSecondaryCode(short *sc, int prn);
void bocMn(const std::vector<short> &tmp_ca, short *ca, int m, int n);

// rinex函数声明
int readContentsData(const std::string &line, double data[39], int data_offset, datetime_t *time);
int readBdsB1CEphemerisCpp(std::vector<ephem_t> eph_vec[MAX_SAT], const std::string &filename);
int matchEpoch(bdstime_t obs_time, vector<ephem_t> eph);
void computeBdtMinMax(const std::vector<ephem_t> eph_vector[MAX_SAT], bdstime_t *bdt_min, bdstime_t *bdt_max, datetime_t *tmin, datetime_t *tmax);
double normalizeBdsAngle(double angle);
int alignBdsAlmanacToa4096(int toa, int *week);
unsigned char getBdsSatTypeFromEphem(const ephem_t *eph);
ephem_t deriveBdsAlmanacFromEphem(const ephem_t *eph, int alm_week, int alm_toa);
void completeBdsAlmanacFromEphem(const ephem_t *eph_list[63], ephem_t alm_out[63], int current_bds_week, int current_bds_tow_seconds);


// gnss_time函数声明
bool isLeapYear(int year);
int daysInMonth(int year, int month);
double subBdsTime(bdstime_t g1, bdstime_t g0);
void dateToBdt(const datetime_t *dt, bdstime_t *bdt);
void bdtToDate(const bdstime_t *bdt, datetime_t *date);
int compareBdt(const bdstime_t *t1, const bdstime_t *t2);
void setScenarioStartTimeBds(bdstime_t *b0, bdstime_t bmin, bdstime_t bmax, datetime_t *t0, datetime_t *tmin, datetime_t *tmax, bool time_overwrite, int neph, vector<ephem_t> eph1[MAX_SAT]);
bdstime_t addBdsTime(bdstime_t g, double dt);
bdstime_t computeSatelliteTxTime(bdstime_t grx, double pseudorange);
b1c_nav_time_fields_t computeB1CNavTimeFields(bdstime_t g);

// bds_sdr函数声明
void *runBdsTask(sim_t *sim);

// channel函数声明
void initChannels(channel_t *chan, vector<int> &allocated_sat);
int allocateChannel(channel_t *chan, vector<ephem_t> *eph_vector, bdstime_t brx, double *xyz, double elv_mask, vector<int> &current_eph_index, vector<int> &allocated_sat, const ephem_t alm[63] = nullptr);

// inav_msg函数声明
int encodeB1CSubframe1Bits(const short nav_bits[14], short encoded_bits[72]);
int generateB1CSubframe1(bdstime_t g, channel_t *chan, ephem_t *eph, short *subframe);
int generateB1CSubframe2(bdstime_t g, channel_t *chan, ephem_t *eph, short *subframe);
int generateB1CSubframe3(bdstime_t g, channel_t *chan, ephem_t *eph, const ephem_t alm[63], short *subframe);
int generateBdsB1CMessage(bdstime_t g, int svid, const ephem_t *eph, const ephem_t alm[63], short nav_msg[1800]);
int generateB1CNavMessage(bdstime_t g, channel_t *chan, ephem_t *eph, const ephem_t alm[63]);

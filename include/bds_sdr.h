/*! \file bds_sdr.h
 *  \brief BDS SDR 主头文件，包含必要的标准库和项目头文件
 *  \author LackWood Du
 *  \date 2025-12-19
 */
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
#include <iomanip>
#include <algorithm>
#include <cstdint>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <map>
#include <complex>
#include <cassert>
#include <time.h>
#include <filesystem>
// #include <uhd/usrp/multi_usrp.hpp>

#include "constants.h"
#include "structures.h"

// geodesy函数声明
bool isGEO(const int &PRN);
void satpos(ephem_t eph, bdstime_t g, double *pos, double *vel, double *clk);
void printSatState(const double *pos, const double *vel, const double *clk);
void subVect(double *y, const double *x1, const double *x2);
double normVect(const double *x);
void xyz2llh(const double *xyz, double *llh);
void llh2xyz(const double *llh, double *xyz);
void ltcmat(const double *llh, double t[3][3]);
void ecef2neu(const double *xyz, double t[3][3], double *neu);
void neu2azel(double *azel, const double *neu);
int checkSatVisibility(ephem_t eph, bdstime_t g, double *xyz, double *azel, int prn, double elvMask);

// bds_sig函数声明
void computeRange(range_t *rho, ephem_t eph, bdstime_t g, double xyz[], int prn);
void computeCodePhase_from_Galileo(channel_t *chan, range_t rho1, double dt, bdstime_t grx);
void computeCodePhase(channel_t *chan, range_t rho1, double dt, bdstime_t grx);
void codegen_B1C_data(short *ca, int prn);
void codegen_B1C_pilot(short *ca11, short *ca61, int prn);
void hex_to_b1c_ca(std::vector<short> &tmp_ca, int prn, int flag);
void BOC11(const std::vector<short> &tmp_ca, short *ca);
void BOC61(const std::vector<short> &tmp_ca, short *ca);
void codegen_B1C_secondary_code(short *sc, int prn);
void BOCmn(const std::vector<short> &tmp_ca, short *ca, int m, int n);
void compareB1C_PRN(int prn, const vector<short> &tmp_ca, const string &csvFile);

// rinex函数声明
void convertD2E(char *line);
int readContentsData(char *str, double *data, datetime_t *time, bool read_time);
int readBdsB1CEphemerisCpp(std::vector<ephem_t> eph_vec[MAX_SAT], const std::string &filename);
void printEphVec(const std::vector<ephem_t> eph_vec[MAX_SAT]);
int epoch_matcher(bdstime_t obsTime, vector<ephem_t> eph);
void compute_bdt_min_max(const std::vector<ephem_t> eph_vector[MAX_SAT], bdstime_t *bdt_min, bdstime_t *bdt_max, datetime_t *tmin, datetime_t *tmax);

// generateCodeFile函数声明
long long modPow(long long base, long long exp, long long mod);
int legendreSymbol(int a, int p);
bool isPrimeInt(int n);
std::vector<int> generateLegendreSequence(int p);
std::vector<int> generateWeilSequence(const std::vector<int> &L, int w);
std::vector<int> generateB1CPrimaryCode(const std::vector<int> &W, int p, int codeLength);
std::string codeToString(const std::vector<int> &code);
void writeOneCodeRow(std::ofstream &file, int prn, const std::vector<int> &code);
std::string chips24ToOctal(const std::vector<int> &chips);
void verifyAgainstICD(int prn, const std::vector<int> &code, const std::string &refHead, const std::string &refTail);
void generateB1CComponentCodes(const std::vector<B1CParams> &params, const std::vector<int> &legendre, const std::string &outputCsvPath, const int codeLength);
std::vector<B1CParams> loadB1CParamsFromCSV(const std::string &csvPath);
void fixOctalFieldWidthInCSV(const std::string &csvPath);

// gnss_time函数声明
bool isLeapYear(int year);
int daysInMonth(int year, int month);
double subBdsTime(bdstime_t g1, bdstime_t g0);
void date2bdt(const datetime_t *dt, bdstime_t *bdt);
void bdt2date(const bdstime_t *bdt, datetime_t *date);
int bdt_cmp(const bdstime_t *t1, const bdstime_t *t2);
void set_scenario_start_time_bds(bdstime_t *b0, bdstime_t bmin, bdstime_t bmax, datetime_t *t0, datetime_t *tmin, datetime_t *tmax, bool timeoverwrite, int neph, vector<ephem_t> eph1[MAX_SAT]);
bdstime_t incBdsTime(bdstime_t g, double dt);

// bds_sdr函数声明
void *bds_task(sim_t *sim);

// channel函数声明
void init_channel(channel_t *chan, vector<int> &allocatedSat);
int allocateChannel(channel_t *chan, vector<ephem_t> *eph_vector, bdstime_t brx, double *xyz, double elvMask, map<int, int> *sm, vector<int> &current_eph_index, vector<int> &allocatedSat);

// inav_msg函数声明
int generateINavMsg(bdstime_t g, channel_t *chan, ephem_t *eph, int idx = 0);
int generateB1CNavMsg(bdstime_t g, channel_t *chan, ephem_t *eph);
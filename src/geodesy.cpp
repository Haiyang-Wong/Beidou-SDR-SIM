/**
 * \file geodesy.cpp
 * \brief GNSS 地球几何计算相关函数实现
 * \author LackWood Du
 * \date 2025-12-19
 */

#include "bds_sdr.h"
#pragma message("Compiling geodesy.cpp")

bool isGEO(const int &PRN)
{
    static const std::vector<int> geoPRNs = {
        01, 02, 03, 04, 05, 59, 60, 61};
    return std::find(geoPRNs.begin(), geoPRNs.end(), PRN) != geoPRNs.end();
}

void satpos(
    ephem_t eph,
    bdstime_t g,
    double *pos,
    double *vel,
    double *clk)
{
    // Computing Satellite Velocity using the Broadcast Ephemeris
    // http://www.ngs.noaa.gov/gps-toolbox/bc_velo.htm

    double tk; // 当前时刻与卫星参考时间的差值
    double mk; // 平近点角
    double ek; // 偏近点角，需要根据kepler方程迭代计算
    double ekold;
    double ekdot;
    double cek, sek;
    double pk; // 真近点角，从椭圆的近地点出发，到卫星当前位置的角度（在轨道平面上测量）
    double pkdot;
    double c2pk, s2pk;
    double uk; // 轨道角，又叫升交距，表示升交点到卫星当前位置的角度（在轨道平面上测量）
    double ukdot;
    double cuk, suk;
    double ok; // 卫星轨道升交点经度
    double sok, cok;
    double ik;
    double ikdot;
    double sik, cik;
    double rk; // 卫星到地球中心的距离
    double rkdot;
    double xpk, ypk;
    double xpkdot, ypkdot;

    double relativistic, OneMinusecosE = 0, tmp;

    tk = g.sec - eph.toe.sec;
    // tk = -14;

    if (tk > SECONDS_IN_HALF_WEEK)
        tk -= SECONDS_IN_WEEK;
    else if (tk < -SECONDS_IN_HALF_WEEK)
        tk += SECONDS_IN_WEEK;

    // cout << "tk = " << tk << endl;

    // 当前时刻的平近点角
    mk = eph.m0 + eph.n * tk;

    // 以下过程通过迭代法，计算偏近点角 ek
    ek = mk;
    ekold = ek + 1.0;

    int cpt = 0;

    // 开普勒方程： M = E - e*sin(E) ——M是平近点角、E是偏近点角、e是轨道离心率
    // 对上述公式对时间T求导之后可得：M' = E' - e*cos(E)*E' ——M'是平均角速度n
    // 因此可以得到： E' = n / (1 - e*cos(E)) ——(1 - e*cos(E))即OneMinusecosE
    while ((fabs(ek - ekold) > 1.0E-14) && cpt < 500)
    {

        cpt++;
        ekold = ek;
        OneMinusecosE = 1.0 - eph.ecc * cos(ekold);
        ek = ek + (mk - ekold + eph.ecc * sin(ekold)) / OneMinusecosE;
    }

    // 计算偏近点角的正弦和余弦值
    sek = sin(ek);
    cek = cos(ek);

    // 计算偏近点角E随时间的变化率
    ekdot = eph.n / OneMinusecosE;

    // 表示相对论改正项，单位是秒，用于修正卫星钟差，ecc是轨道离心率，sqrta是轨道半长轴的平方根，sek是偏近点角的正弦值
    relativistic = -4.442807633E-10 * eph.ecc * eph.sqrta * sek;

    // 轨道角aop是近地点幅角，是升交点到近地点的角度（轨道平面上计算）
    pk = atan2(eph.sq1e2 * sek, cek - eph.ecc) + eph.aop;
    pkdot = eph.sq1e2 * ekdot / OneMinusecosE;

    s2pk = sin(2.0 * pk);
    c2pk = cos(2.0 * pk);

    // 计算卫星的修正轨道角
    uk = pk + eph.cus * s2pk + eph.cuc * c2pk;
    // 轨道角的正余弦值
    suk = sin(uk);
    cuk = cos(uk);
    // 轨道角的变化率
    ukdot = pkdot * (1.0 + 2.0 * (eph.cus * c2pk - eph.cuc * s2pk));

    // 计算修正轨道半径及其变化率
    rk = eph.A * OneMinusecosE + eph.crc * c2pk + eph.crs * s2pk;
    rkdot = eph.A * eph.ecc * sek * ekdot + 2.0 * pkdot * (eph.crs * c2pk - eph.crc * s2pk);

    // 计算卫星轨道倾角和其变化率
    ik = eph.inc0 + eph.idot * tk + eph.cic * c2pk + eph.cis * s2pk;
    sik = sin(ik);
    cik = cos(ik);
    ikdot = eph.idot + 2.0 * pkdot * (eph.cis * c2pk - eph.cic * s2pk);

    // 计算在卫星轨道平面内，卫星的坐标及其速度分量
    xpk = rk * cuk;
    ypk = rk * suk;
    xpkdot = rkdot * cuk - ypk * ukdot;
    ypkdot = rkdot * suk + xpk * ukdot;

    // 计算升交点经度及其正余弦值
    ok = eph.omg0 + tk * eph.omgkdot - OMEGA_EARTH * eph.toe.sec;
    sok = sin(ok);
    cok = cos(ok);

    // Compute position
    // 卫星位置，ECEF坐标系中的x、y、z坐标
    pos[0] = xpk * cok - ypk * cik * sok;
    pos[1] = xpk * sok + ypk * cik * cok;
    pos[2] = ypk * sik;

    tmp = ypkdot * cik - ypk * sik * ikdot;

    // Compute velocity
    // 卫星速度，ECEF坐标系中的vx、vy、vz分量
    vel[0] = -eph.omgkdot * pos[1] + xpkdot * cok - tmp * sok;
    vel[1] = eph.omgkdot * pos[0] + xpkdot * sok + tmp * cok;
    vel[2] = ypk * cik * ikdot + ypkdot * sik;

    if (isGEO(eph.PRN))
    { // 请根据您的实际数据结构修改此判断条件
        // cout << "Applying GEO pole tide correction for PRN " << eph.PRN << endl;
        double fi = OMEGA_EARTH * tk;     // 计算地球自转角度修正参数
        double five = 180.0 / K_PI * 5.0; // 5度偏差角，转换为弧度

        // 计算旋转所需的正余弦值
        double cos_fi = cos(fi);
        double sin_fi = sin(fi);
        double cos_five = cos(five);
        double sin_five = sin(five);

        // 提取改正前的坐标
        double X_old = pos[0];
        double Y_old = pos[1];
        double Z_old = pos[2];

        // 应用极移改正矩阵（旋转顺序：先绕Z轴旋转 -five，再绕X轴旋转 fi）
        // 这对应于 Python 代码中的 rx(fi) * rz(-five) * a.T
        pos[0] = cos_five * X_old - sin_five * Y_old;
        pos[1] = sin_five * cos_fi * X_old + cos_five * cos_fi * Y_old + sin_fi * Z_old;
        pos[2] = -sin_five * sin_fi * X_old - cos_five * sin_fi * Y_old + cos_fi * Z_old;

        // 注意：如果速度信息也需要同样的改正，请对 vel[0], vel[1], vel[2] 进行相同的矩阵变换
        double vx_old = vel[0];
        double vy_old = vel[1];
        double vz_old = vel[2];

        vel[0] = cos_five * vx_old - sin_five * vy_old;
        vel[1] = sin_five * cos_fi * vx_old + cos_five * cos_fi * vy_old + sin_fi * vz_old;
        vel[2] = -sin_five * sin_fi * vx_old - cos_five * sin_fi * vy_old + cos_fi * vz_old;
    }

    // Satellite clock correction
    // 计算卫星钟差
    tk = g.sec - eph.toc.sec;

    if (tk > SECONDS_IN_HALF_WEEK)
        tk -= SECONDS_IN_WEEK;
    else if (tk < -SECONDS_IN_HALF_WEEK)
        tk += SECONDS_IN_WEEK;

    // 卫星的总钟差，即卫星钟多项式修正 + 相对论改正项 relativistic - 广播星历给定的 bgde5b 偏置
    // af0：常数项偏差，表示卫星在参考时间 toc 的零时刻钟差
    // af1 * tk：一阶项，卫星钟漂移率
    // af2 * tk^2：二阶项，卫星钟加速度
    // relativistic：相对论效应引起的钟差修正
    // bgde5b：广播星历中给定的偏置值
    // 卫星时钟相较于系统时钟的误差
    clk[0] = eph.af0 + tk * (eph.af1 + tk * eph.af2) + relativistic;

    // 卫星钟差速率，上述公式对时间 tk 求导得到
    clk[1] = eph.af1 + 2.0 * tk * eph.af2;

    return;
}

void printSatState(const double *pos, const double *vel, const double *clk)
{
    std::cout << std::fixed << std::setprecision(6);

    std::cout << "================ Satellite State ================\n";

    std::cout << "Position (ECEF) [m]\n";
    std::cout << "  X = " << std::setw(15) << pos[0] << "\n";
    std::cout << "  Y = " << std::setw(15) << pos[1] << "\n";
    std::cout << "  Z = " << std::setw(15) << pos[2] << "\n\n";

    std::cout << "Velocity (ECEF) [m/s]\n";
    std::cout << "  VX = " << std::setw(15) << vel[0] << "\n";
    std::cout << "  VY = " << std::setw(15) << vel[1] << "\n";
    std::cout << "  VZ = " << std::setw(15) << vel[2] << "\n\n";

    std::cout << "Clock\n";
    std::cout << "  Clock Bias   (s)  = " << std::setw(15) << clk[0] << "\n";
    std::cout << "  Clock Drift  (s/s)= " << std::setw(15) << clk[1] << "\n";

    std::cout << "=================================================\n";
}

/**
 * @brief 向量减法
 */
void subVect(double *y, const double *x1, const double *x2)
{
    y[0] = x1[0] - x2[0];
    y[1] = x1[1] - x2[1];
    y[2] = x1[2] - x2[2];

    return;
}

/**
 * @brief 计算向量模长
 */
double normVect(const double *x)
{
    return (sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]));
}

/**
 * @brief 主要功能是将地心地固坐标系中的坐标（X，Y，Z）转为经纬度高度（Latitude，Longitude，Height）
 */
void xyz2llh(const double *xyz, double *llh)
{
    double a, eps, e, e2;
    double x, y, z;
    double rho2, dz, zdz, nh, slat, n, dz_new;

    a = WGS84_RADIUS;
    e = WGS84_ECCENTRICITY;

    eps = 1.0e-3;
    e2 = e * e;

    if (normVect(xyz) < eps)
    {
        // Invalid ECEF vector
        llh[0] = 0.0;
        llh[1] = 0.0;
        llh[2] = -a;

        return;
    }

    x = xyz[0];
    y = xyz[1];
    z = xyz[2];

    rho2 = x * x + y * y;
    dz = e2 * z;

    while (1)
    {
        zdz = z + dz;
        nh = sqrt(rho2 + zdz * zdz);
        slat = zdz / nh;
        n = a / sqrt(1.0 - e2 * slat * slat);
        dz_new = n * e2 * slat;

        if (fabs(dz - dz_new) < eps)
            break;

        dz = dz_new;
    }

    llh[0] = atan2(zdz, sqrt(rho2));
    llh[1] = atan2(y, x);
    llh[2] = nh - n;

    return;
}

/**
 * \brief 主要功能是将经纬度高度（Latitude，Longitude，Height）转为地心地固坐标系中的坐标（X，Y，Z）
 * \param[in] llh 长度为3的数组，顺序是{纬度、经度、高度}，单位为弧度和米
 * \param[out] xyz 长度为3的数组，顺序是{X、Y、Z}，单位为米
 */
void llh2xyz(const double *llh, double *xyz)
{
    double n;
    double a;
    double e;
    double e2;
    double clat;
    double slat;
    double clon;
    double slon;
    double d, nph;
    double tmp;

    a = WGS84_RADIUS;
    e = WGS84_ECCENTRICITY;
    e2 = e * e;

    clat = cos(llh[0]);
    slat = sin(llh[0]);
    clon = cos(llh[1]);
    slon = sin(llh[1]);
    d = e * slat;

    n = a / sqrt(1.0 - d * d);
    nph = n + llh[2];

    tmp = nph * clat;
    xyz[0] = tmp * clon;
    xyz[1] = tmp * slon;
    xyz[2] = ((1.0 - e2) * n + llh[2]) * slat;

    return;
}

/**
 * @brief 计算局部切平面转换矩阵
 * @note 局部切平面转换矩阵可用于将ECEF坐标系转换为局部坐标系（如东北天坐标系NEU）
 */
void ltcmat(const double *llh, double t[3][3])
{
    double slat, clat;
    double slon, clon;

    slat = sin(llh[0]);
    clat = cos(llh[0]);
    slon = sin(llh[1]);
    clon = cos(llh[1]);

    t[0][0] = -slat * clon;
    t[0][1] = -slat * slon;
    t[0][2] = clat;
    t[1][0] = -slon;
    t[1][1] = clon;
    t[1][2] = 0.0;
    t[2][0] = clat * clon;
    t[2][1] = clat * slon;
    t[2][2] = slat;

    return;
}

/**
 * @brief 将地心地固坐标系转为东北天坐标系
 * \param[in] xyz 地心地固坐标系下的坐标
 * \param[in] t 局部切平面转换矩阵
 * \param[out] neu 东北天坐标系下的坐标
 */
void ecef2neu(const double *xyz, double t[3][3], double *neu)
{
    neu[0] = t[0][0] * xyz[0] + t[0][1] * xyz[1] + t[0][2] * xyz[2];
    neu[1] = t[1][0] * xyz[0] + t[1][1] * xyz[1] + t[1][2] * xyz[2];
    neu[2] = t[2][0] * xyz[0] + t[2][1] * xyz[1] + t[2][2] * xyz[2];

    return;
}

/**
 * @brief 将东北天坐标系转为方位角和俯仰角
 */
void neu2azel(double *azel, const double *neu)
{
    double ne;

    azel[0] = atan2(neu[1], neu[0]);
    if (azel[0] < 0.0)
        azel[0] += (2.0 * K_PI);

    ne = sqrt(neu[0] * neu[0] + neu[1] * neu[1]);
    azel[1] = atan2(neu[2], ne);

    return;
}

/**
 * @brief 计算卫星的方位角和俯仰角，并于给定阈值比较，检查卫星是否可见
 * \param[in] eph 卫星星历数据
 * \param[in] g 当前时间
 * \param[in] xyz 接收机位置坐标
 * \param[in] elvMask（Elevation angle），单位是度，卫星与接收机连线处于水平时为0度，垂直向上时为90度，该变量用于屏蔽低仰角卫星
 * \param[out] azel 卫星方位角和仰角
 * \param[in] prn 卫星编号
 * \note 计算卫星相较于接收机的方位角和俯仰角，并根据设定的仰角掩膜判断卫星是否可见
 * \return 返回值：1表示卫星可见，0表示卫星不可见，-1表示星历数据无效
 */
int checkSatVisibility(
    ephem_t eph, bdstime_t g,
    double *xyz,
    double *azel,
    int prn,
    double elvMask = 10) // checked
{
    double llh[3], neu[3];
    double pos[3], vel[3], clk[3], los[3];
    double tmat[3][3];
    // 判断星历数据是否有效
    if (eph.vflg != 1)
    {
        return (-1); // Invalid}
    }
    xyz2llh(xyz, llh);
    ltcmat(llh, tmat);
    // 计算卫星位置、速度和时钟偏差，返回值分别存储在 pos、vel 和 clk 数组中
    satpos(eph, g, pos, vel, clk);
    // los = pos - xyz
    subVect(los, pos, xyz);
    // 将ECEF中的los，根据tmat转换为NEU坐标系（东北天坐标系）
    ecef2neu(los, tmat, neu);
    // 根据NEU坐标系下的los，计算卫星的方位角和仰角
    neu2azel(azel, neu);
    // 判断卫星是否可见
    if (azel[1] * R2D > elvMask)
        return (1); // Visible
    // else
    return (0); // Invisible
}
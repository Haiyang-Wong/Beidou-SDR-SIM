#include "bds_sdr.h"

bool isGEO(const int &prn)
{
    static const std::vector<int> geo_prns = {
        01, 02, 03, 04, 05, 59, 60, 61};
    return std::find(geo_prns.begin(), geo_prns.end(), prn) != geo_prns.end();
}

void satPos(
    ephem_t eph,
    bdstime_t g,
    double *pos,
    double *vel,
    double *clk)
{

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

    double relativistic, one_minus_ecos_e = 0, tmp;

    const bdstime_t toe_time = bdsTimeFromWeekSeconds(eph.week, eph.toe);
    tk = subBdsTime(g, toe_time);
    if (tk > SECONDS_IN_HALF_WEEK)
        tk -= SECONDS_IN_WEEK;
    else if (tk < -SECONDS_IN_HALF_WEEK)
        tk += SECONDS_IN_WEEK;
    mk = eph.m0 + eph.n * tk;
    ek = mk;
    ekold = ek + 1.0;

    int cpt = 0;
    while ((fabs(ek - ekold) > 1.0E-14) && cpt < 500)
    {

        cpt++;
        ekold = ek;
        one_minus_ecos_e = 1.0 - eph.ecc * cos(ekold);
        ek = ek + (mk - ekold + eph.ecc * sin(ekold)) / one_minus_ecos_e;
    }
    sek = sin(ek);
    cek = cos(ek);
    ekdot = eph.n / one_minus_ecos_e;
    relativistic = -4.442807633E-10 * eph.ecc * eph.sqrt_a * sek;
    pk = atan2(eph.root_ecc * sek, cek - eph.ecc) + eph.w;
    pkdot = eph.root_ecc * ekdot / one_minus_ecos_e;

    s2pk = sin(2.0 * pk);
    c2pk = cos(2.0 * pk);
    uk = pk + eph.cus * s2pk + eph.cuc * c2pk;
    suk = sin(uk);
    cuk = cos(uk);
    ukdot = pkdot * (1.0 + 2.0 * (eph.cus * c2pk - eph.cuc * s2pk));
    rk = eph.axis * one_minus_ecos_e + eph.crc * c2pk + eph.crs * s2pk;
    rkdot = eph.axis * eph.ecc * sek * ekdot + 2.0 * pkdot * (eph.crs * c2pk - eph.crc * s2pk);
    ik = eph.i0 + eph.idot * tk + eph.cic * c2pk + eph.cis * s2pk;
    sik = sin(ik);
    cik = cos(ik);
    ikdot = eph.idot + 2.0 * pkdot * (eph.cis * c2pk - eph.cic * s2pk);
    xpk = rk * cuk;
    ypk = rk * suk;
    xpkdot = rkdot * cuk - ypk * ukdot;
    ypkdot = rkdot * suk + xpk * ukdot;
    ok = eph.omega0 + tk * eph.omega_delta - OMEGA_EARTH * eph.toe;
    sok = sin(ok);
    cok = cos(ok);
    pos[0] = xpk * cok - ypk * cik * sok;
    pos[1] = xpk * sok + ypk * cik * cok;
    pos[2] = ypk * sik;

    tmp = ypkdot * cik - ypk * sik * ikdot;
    vel[0] = -eph.omega_delta * pos[1] + xpkdot * cok - tmp * sok;
    vel[1] = eph.omega_delta * pos[0] + xpkdot * sok + tmp * cok;
    vel[2] = ypk * cik * ikdot + ypkdot * sik;

    if (isGEO(eph.svid))
    {
        double fi = OMEGA_EARTH * tk;     // 计算地球自转角度修正参数
        double five = 180.0 / K_PI * 5.0; // 5度偏差角，转换为弧度
        double cos_fi = cos(fi);
        double sin_fi = sin(fi);
        double cos_five = cos(five);
        double sin_five = sin(five);
        double x_old = pos[0];
        double y_old = pos[1];
        double z_old = pos[2];
        pos[0] = cos_five * x_old - sin_five * y_old;
        pos[1] = sin_five * cos_fi * x_old + cos_five * cos_fi * y_old + sin_fi * z_old;
        pos[2] = -sin_five * sin_fi * x_old - cos_five * sin_fi * y_old + cos_fi * z_old;
        double vx_old = vel[0];
        double vy_old = vel[1];
        double vz_old = vel[2];

        vel[0] = cos_five * vx_old - sin_five * vy_old;
        vel[1] = sin_five * cos_fi * vx_old + cos_five * cos_fi * vy_old + sin_fi * vz_old;
        vel[2] = -sin_five * sin_fi * vx_old - cos_five * sin_fi * vy_old + cos_fi * vz_old;
    }
    const bdstime_t toc_time = bdsTimeFromWeekSeconds(eph.week, eph.toc);
    tk = subBdsTime(g, toc_time);

    if (tk > SECONDS_IN_HALF_WEEK)
        tk -= SECONDS_IN_WEEK;
    else if (tk < -SECONDS_IN_HALF_WEEK)
        tk += SECONDS_IN_WEEK;
    clk[0] = eph.af0 + tk * (eph.af1 + tk * eph.af2) + relativistic;
    clk[1] = eph.af1 + 2.0 * tk * eph.af2;

    return;
}


/**
 * @brief 向量减法
 */
void subVector(double *y, const double *x1, const double *x2)
{
    y[0] = x1[0] - x2[0];
    y[1] = x1[1] - x2[1];
    y[2] = x1[2] - x2[2];

    return;
}

/**
 * @brief 计算向量模长
 */
double normVector(const double *x)
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

    if (normVector(xyz) < eps)
    {
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
 * \param[in] elv_mask（Elevation angle），单位是度，卫星与接收机连线处于水平时为0度，垂直向上时为90度，该变量用于屏蔽低仰角卫星
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
    double elv_mask = 10)
{
    double llh[3], neu[3];
    double pos[3], vel[3], clk[3], los[3];
    double tmat[3][3];
    if (eph.valid != 1)
    {
        return (-1); // Invalid}
    }
    xyz2llh(xyz, llh);
    ltcmat(llh, tmat);
    satPos(eph, g, pos, vel, clk);
    subVector(los, pos, xyz);
    ecef2neu(los, tmat, neu);
    neu2azel(azel, neu);
    if (azel[1] * R2D > elv_mask)
        return (1); // Visible
    return (0); // Invisible
}

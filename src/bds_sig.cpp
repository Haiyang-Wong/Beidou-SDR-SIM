#include "bds_sdr.h"

/*! \brief 计算卫星和接收机之间的伪距
 *  \param[out] rho 其中包含伪距、方位角、仰角等信息的结构体
 *  \param[in] eph 卫星的星历数据
 *  \param[in] g 当前的Galileo时间
 *  \param[in] xyz 接收机的位置
 */
void computeRange(
    range_t *rho,
    ephem_t eph,
    bdstime_t g,
    double xyz[],
    int prn)
{
    double pos[3] = {0.0}, vel[3] = {0.0}, clk[2] = {0.0};
    double los[3] = {0.0};
    double tau = 0.0;
    double range = 0.0;
    double xrot = 0.0, yrot = 0.0;

    double llh[3] = {0.0}, neu[3] = {0.0};
    double tmat[3][3] = {0};
    satPos(eph, g, pos, vel, clk);
    subVector(los, pos, xyz);
    tau = normVector(los) / SPEED_OF_LIGHT;
    pos[0] -= vel[0] * tau;
    pos[1] -= vel[1] * tau;
    pos[2] -= vel[2] * tau;
    xrot = pos[0] + pos[1] * GNSS_OMEGA_EARTH_DOT * tau;
    yrot = pos[1] - pos[0] * GNSS_OMEGA_EARTH_DOT * tau;
    pos[0] = xrot;
    pos[1] = yrot;
    subVector(los, pos, xyz);
    range = normVector(los);
    rho->d = range;
    rho->range = range - SPEED_OF_LIGHT * clk[0];
    xyz2llh(xyz, llh);    // convert userXYZ to llh
    ltcmat(llh, tmat);
    ecef2neu(los, tmat, neu);
    neu2azel(rho->azel, neu);
    rho->g = g;

    return;
}

/*! \brief 计算指定卫星通道的码相位
 *  \param[in,out] chan 正在处理的通道指针（函数内会更新其状态）
 *  \param[in] rho1 在经过 dt 时间后的当前伪距等量测
 *  \param[in] dt 时间增量（单位：秒）
 *  \param[in] grx 接收机的Galileo时间
 */
void computeCodePhase(
    channel_t *chan,
    range_t rho1,
    double dt,
    bdstime_t grx)
{
    double rhorate;
    rhorate = (rho1.range - chan->rho0.range) / dt;
    chan->f_carr = (-rhorate / LAMBDA_B1C);

    chan->f_code = CODE_FREQ_B1C + chan->f_carr * CARR_TO_CODE_B1C;
    if (chan->set_code_phase)
    {
        bdstime_t gtx = computeSatelliteTxTime(grx, rho1.range);
        constexpr double eps = 1e-6;
        double ms = (static_cast<double>(gtx.milliseconds) + gtx.sub_milliseconds);

        chan->code_phase = (fmod(ms, 10.0) / 10.0) * CA_SEQ_LEN_B1C;

        const b1c_nav_time_fields_t time_fields = computeB1CNavTimeFields(gtx);
        double frame_offset = ((static_cast<double>(gtx.milliseconds) + gtx.sub_milliseconds) - time_fields.frame_start_ms) / BDS_MILLISECONDS_IN_SECOND;
        if (frame_offset < 0.0 && fabs(frame_offset) < eps)
            frame_offset = 0.0;
        frame_offset = fmod(frame_offset + eps, 18.0);
        if (frame_offset < 0.0)
            frame_offset += 18.0;

        chan->ibit = static_cast<int>(floor(frame_offset * 100.0));
        if (chan->ibit >= NAV_SEQ_LEN_B1C)
            chan->ibit = 0;
        chan->set_code_phase = false;
    }
    chan->rho0 = rho1;
    return;
}

/**
 * \brief 将十六进制字符串转换为B1C码片序列
 * \param[out] tmp_ca 存储转换后B1C码片序列
 * \param[in] prn 卫星PRN号
 * \param[in] flag 指定转换的码类型（主码或辅码），取值为 B1C_DATA_PRIMARY、B1C_PILOT_PRIMARY、B1C_PILOT_SUB
 */
void hexToB1CCa(
    std::vector<short> &tmp_ca,
    int prn,
    int flag)
{
    const std::string *hex_str = nullptr;
    int target_len = 0;
    switch (flag)
    {
    case B1C_DATA_PRIMARY:
        hex_str = &b1c_data_primary_hex_ca.at(prn - 1);
        target_len = CA_SEQ_LEN_B1C;
        break;

    case B1C_PILOT_PRIMARY:
        hex_str = &b1c_pilot_primary_hex_ca.at(prn - 1);
        target_len = CA_SEQ_LEN_B1C;
        break;

    case B1C_PILOT_SUB:
        hex_str = &b1c_pilot_sub_hex_ca.at(prn - 1);
        target_len = SC_SEQ_LEN_B1C;
        break;

    default:
        throw std::invalid_argument("Invalid B1C flag");
    }

    tmp_ca.clear();
    tmp_ca.reserve(target_len);
    for (char ch : *hex_str)
    {
        if (tmp_ca.size() >= static_cast<size_t>(target_len))
        {
            break;
        }

        char c = std::toupper(static_cast<unsigned char>(ch));
        int value = 0;

        if (c >= '0' && c <= '9')
        {
            value = c - '0';
        }
        else if (c >= 'A' && c <= 'F')
        {
            value = c - 'A' + 10;
        }
        else
        {
            continue;
        }
        for (int i = 3; i >= 0; --i)
        {
            if (tmp_ca.size() >= static_cast<size_t>(target_len))
            {
                break;
            }

            int bit = (value >> i) & 0x1;
            tmp_ca.push_back(bit ? -1 : 1);
        }
    }
}


void bocMn(
    const std::vector<short> &tmp_ca,
    short *ca,
    int m, // 副载波倍数
    int n  // 码速率倍数（B1C中 n=1）
)
{
    const int chips_per_code = tmp_ca.size();
    const int subchips_per_chip = 2 * m; // 半周期数
    int sc = 0;                          // 全局子码片计数器

    for (int chip = 0; chip < chips_per_code; ++chip)
    {
        for (int k = 0; k < subchips_per_chip; ++k)
        {
            int subcarrier = ((sc & 1) == 0) ? -1 : 1;
            ca[chip * subchips_per_chip + k] = tmp_ca[chip] * subcarrier;
            sc++;
        }
    }
}


void generateB1CDataCode(short *ca, int prn)
{
    vector<short> tmp_ca;
    hexToB1CCa(tmp_ca, prn, B1C_DATA_PRIMARY);
    bocMn(tmp_ca, ca, 1, 1);
}

void generateB1CPilotCode(short *ca11, short *ca61, int prn)
{
    vector<short> tmp_ca;
    hexToB1CCa(tmp_ca, prn, B1C_PILOT_PRIMARY);
    bocMn(tmp_ca, ca11, 1, 1);
    bocMn(tmp_ca, ca61, 6, 1);
}

void generateB1CSecondaryCode(short *sc, int prn)
{
    vector<short> tmp_ca;
    hexToB1CCa(tmp_ca, prn, B1C_PILOT_SUB);
    for (size_t i = 0; i < tmp_ca.size(); ++i)
    {
        sc[i] = tmp_ca[i];
    }
}

#include "bds_sdr.h"

namespace
{
void updateWeightedB1CCode(channel_t &chan)
{
    for (int i = 0; i < 2 * CA_SEQ_LEN_B1C; ++i)
    {
        chan.weighted_b1c_data[i] = B1C_DATA_WEIGHT * chan.ca_b1c_data[i];
        chan.weighted_b1c_pilot11[i] = B1C_PILOT11_WEIGHT * chan.ca_b1c_pilot11[i];
    }

    for (int i = 0; i < 12 * CA_SEQ_LEN_B1C; ++i)
        chan.weighted_b1c_pilot61[i] = B1C_PILOT61_WEIGHT * chan.ca_b1c_pilot61[i];
}
}

/**
 * @brief 初始化通道数组、清空卫星分配状态, 主要完成动态内存分配和状态初始化工作
 * @param chan 指向通道数组的指针，通常用于表示多个卫星信号的生成通道
 * @param allocated_sat 指向整型数组的指针，用于某个卫星被分配到了哪一个通道，或者是哪个卫星对应哪个通道
 */
void initChannels(channel_t *chan, vector<int> &allocated_sat)
{
    for (int i = 0; i < MAX_CHAN; i++)
    {
        chan[i].set_code_phase = true;                                                 // 代码相位设置标志
        chan[i].prn = 0;                                                               // 将PRN号初始化为0
        chan[i].ca_b1c_data = (short *)malloc(2 * CA_SEQ_LEN_B1C * sizeof(short));     // 为B1C数据分量分配内存
        chan[i].ca_b1c_pilot11 = (short *)malloc(2 * CA_SEQ_LEN_B1C * sizeof(short));  // 为B1C导频分量分配内存
        chan[i].ca_b1c_pilot61 = (short *)malloc(12 * CA_SEQ_LEN_B1C * sizeof(short));
        chan[i].weighted_b1c_data = (double *)malloc(2 * CA_SEQ_LEN_B1C * sizeof(double));
        chan[i].weighted_b1c_pilot11 = (double *)malloc(2 * CA_SEQ_LEN_B1C * sizeof(double));
        chan[i].weighted_b1c_pilot61 = (double *)malloc(12 * CA_SEQ_LEN_B1C * sizeof(double)); // 为B1C导频分量分配内存
        chan[i].b1c_pilot_sub_code = (short *)malloc(SC_SEQ_LEN_B1C * sizeof(short));  // 为次级码分配内存
        chan[i].nav_bit = (short *)malloc(NAV_SEQ_LEN_B1C * sizeof(short));            // 为导航bit分配内存
    }
    for (int sv = 0; sv < MAX_SAT; sv++)
        allocated_sat[sv] = -1;
}

/**
 * \brief 为卫星分配通道
 */
int allocateChannel(channel_t *chan,
                    vector<ephem_t> *eph_vector, // 卫星的星历数据
                    bdstime_t brx,               // 当前仿真时间，以Galileo时间表示
                    double *xyz,                 // ECEF坐标系下，接收机的位置坐标
                    double elv_mask,              //
                    vector<int> &current_eph_index, // 当前每颗卫星使用的星历索引
                    vector<int> &allocated_sat,
                    const ephem_t alm[63])      // 记录每个卫星被分配到了哪个通道
{
    int nsat = 0;
    int i, sv;
    double azel[2];

    range_t rho;
    range_t ref_rho;
    double ref[3] = {0.0};
    double r_ref, r_xyz;
    double phase_ini;
    bool printed_initial_header = false;

    ephem_t eph;
    for (sv = 0; sv < MAX_SAT; sv++)
    {
        if (!eph_vector[sv].size())
            continue;
        if (!eph_vector[sv][0].valid)
            continue;
        current_eph_index[sv] = matchEpoch(brx, eph_vector[sv]);
        if (current_eph_index[sv] < 0)
            continue;
        eph = eph_vector[sv][current_eph_index[sv]];
        if (checkSatVisibility(eph, brx, xyz, azel, sv + 1, elv_mask) == 1)
        {

            nsat++; // 可见卫星的数量加1
            if (allocated_sat[sv] == -1) // Visible but not allocated
            {
                for (i = 0; i < MAX_CHAN; i++)
                {
                    if (chan[i].prn == 0)
                    {
                        chan[i].prn = sv + 1;
                        chan[i].azel[0] = azel[0];
                        chan[i].azel[1] = azel[1];
                        chan[i].set_code_phase = true; // 要求下次 computeCodePhase 校准码相位和 ibit
                        generateB1CDataCode(chan[i].ca_b1c_data, chan[i].prn);
                        generateB1CPilotCode(chan[i].ca_b1c_pilot11, chan[i].ca_b1c_pilot61, chan[i].prn);
                        updateWeightedB1CCode(chan[i]);
                        generateB1CSecondaryCode(chan[i].b1c_pilot_sub_code, chan[i].prn);
                        computeRange(&rho, eph, brx, xyz, chan[i].prn);
                        chan[i].rho0 = rho;
                        r_xyz = rho.range;

                        computeRange(&ref_rho, eph, brx, ref, chan[i].prn);
                        r_ref = ref_rho.range;
                        phase_ini = (2.0 * r_ref - r_xyz) / LAMBDA_B1C;
                        chan[i].carr_phase = phase_ini - floor(phase_ini);
                        bdstime_t g_tx = computeSatelliteTxTime(brx, chan[i].rho0.range);
                        double ms = (static_cast<double>(g_tx.milliseconds) + g_tx.sub_milliseconds);
                        if (fmod(ms, 10.0) < 0.0)
                            ms += 10.0;
                        double initial_code_phase = (fmod(ms, 10.0) / 10.0) * CA_SEQ_LEN_B1C;
                        generateB1CNavMessage(g_tx, &chan[i], &eph, alm);
                        if (!printed_initial_header)
                        {
                            cout << "------------------ Initial visible satellite signal states -------------------" << endl;
                            printed_initial_header = true;
                        }
                        cout << "  prn" << chan[i].prn
                             << "  tx_time=" << g_tx.week << ":" << fixed << setprecision(6) << ((static_cast<double>(g_tx.milliseconds) + g_tx.sub_milliseconds) / BDS_MILLISECONDS_IN_SECOND)
                             << "  code_phase=" << setprecision(3) << initial_code_phase << " chips"
                             << "  pseudorange=" << setprecision(3) << chan[i].rho0.range << " m"
                             << defaultfloat << setprecision(6) << endl;
                        break;
                    }
                }
                if (i < MAX_CHAN)
                    allocated_sat[sv] = i;
            }
        }
        else if (allocated_sat[sv] >= 0) // Not visible but allocated
        {
            chan[allocated_sat[sv]].prn = 0;
            allocated_sat[sv] = -1;
        }
    }
    return (nsat);
}

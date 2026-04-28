/**
 * @file channel.cpp
 * @brief GNSS 信道处理相关函数实现
 * @author LackWood Du
 * @date 2025-12-24
 */
#include "bds_sdr.h"

/**
 * @brief 初始化通道数组、清空卫星分配状态, 主要完成动态内存分配和状态初始化工作
 * @param chan 指向通道数组的指针，通常用于表示多个卫星信号的生成通道
 * @param allocatedSat 指向整型数组的指针，用于某个卫星被分配到了哪一个通道，或者是哪个卫星对应哪个通道
 */
void init_channel(channel_t *chan, vector<int> &allocatedSat)
{
    // Clear all channels
    for (int i = 0; i < MAX_CHAN; i++)
    {
        // Set length of code array at runtime based on samples
        chan[i].set_code_phase = true;                                                 // 代码相位设置标志
        chan[i].prn = 0;                                                               // 将PRN号初始化为0
        chan[i].ca_B1C_data = (short *)malloc(2 * CA_SEQ_LEN_B1C * sizeof(short));     // 为B1C数据分量分配内存
        chan[i].ca_B1C_pilot11 = (short *)malloc(2 * CA_SEQ_LEN_B1C * sizeof(short));  // 为B1C导频分量分配内存
        chan[i].ca_B1C_pilot61 = (short *)malloc(12 * CA_SEQ_LEN_B1C * sizeof(short)); // 为B1C导频分量分配内存
        chan[i].B1C_pilot_sub_code = (short *)malloc(SC_SEQ_LEN_B1C * sizeof(short));  // 为次级码分配内存
        chan[i].nav_bit = (short *)malloc(NAV_SEQ_LEN_B1C * sizeof(short));            // 为导航bit分配内存
        // chan[i].page = (int *)malloc(PAGE_SIZE * sizeof(int)); // 为页面数据分配内存，B1C信号中应该没有页面这个结构
    }

    // Clear satellite allocation flag
    // 卫星分配状态标记为-1, 表示该卫星尚未被分配给任何通道，数组索引表示卫星编号，值表示分配的通道号
    for (int sv = 0; sv < MAX_SAT; sv++)
        allocatedSat[sv] = -1;
}

/**
 * \brief 为卫星分配通道
 */
int allocateChannel(channel_t *chan,
                    // map<int, vector<Rinex3Nav::DataGAL>> *navGAL,
                    vector<ephem_t> *eph_vector, // 卫星的星历数据
                    bdstime_t brx,               // 当前仿真时间，以Galileo时间表示
                    double *xyz,                 // ECEF坐标系下，接收机的位置坐标
                    double elvMask,              //
                    map<int, int> *sm,
                    vector<int> &current_eph_index, // 当前每颗卫星使用的星历索引
                    vector<int> &allocatedSat)      // 记录每个卫星被分配到了哪个通道
{
    // cout << "Allocating Channels..." << endl;
    int nsat = 0;
    int i, sv;
    double azel[2];

    range_t rho;
    double ref[3] = {0.0};
    double r_ref, r_xyz;
    double phase_ini;

    ephem_t eph;

    // 遍历每颗卫星
    for (sv = 0; sv < MAX_SAT; sv++)
    {
        // 判断当前卫星是否有星历，如果没有则跳过
        if (!eph_vector[sv].size())
            continue;

        // 当前卫星的星历是否有效，如果无效则跳过
        if (!eph_vector[sv][0].vflg)
            continue;

        datetime_t t;
        bdt2date(&brx, &t);

        // 找到合适的星历数据，并将其索引保存到 current_eph_index[sv]
        current_eph_index[sv] = epoch_matcher(brx, eph_vector[sv]);
        // cout << "here : " << current_eph_index[sv] << " : " << sv<<  endl;
        // 如果索引值为负数，表示没有找到合适的星历
        if (current_eph_index[sv] < 0)
            continue;

        // 获取星历数据
        eph = eph_vector[sv][current_eph_index[sv]];
        // cout << "GRX: " << grx.sec << " : " << current_eph_index[sv] << endl;

        // 判断卫星是否可见，返回值为1表示可见，0表示不可见
        if (checkSatVisibility(eph, brx, xyz, azel, sv + 1, 10) == 1)
        {

            nsat++; // 可见卫星的数量加1

            // 当前可见的卫星还没有被分配通道
            if (allocatedSat[sv] == -1) // Visible but not allocated
            {
                // Allocated new satellite
                // 遍历所有通道，为上述可见卫星分配一个空闲通道
                for (i = 0; i < MAX_CHAN; i++)
                {
                    // 当前通道尚未被分配卫星
                    if (chan[i].prn == 0)
                    {

                        // Initialize channel
                        chan[i].prn = sv + 1;
                        chan[i].azel[0] = azel[0];
                        chan[i].azel[1] = azel[1];
                        chan[i].g0 = brx; // 分配通道的时间是卫星信号开始观测的时间

                        // Insert latest channel assignment to the map
                        // sm->insert({chan[i].prn, i});

                        // C/A code generation
                        // 根据PRN号，查找固定的主码，并进行SBOC调制，得到最终的B1C码片序列
                        codegen_B1C_data(chan[i].ca_B1C_data, chan[i].prn);
                        codegen_B1C_pilot(chan[i].ca_B1C_pilot11, chan[i].ca_B1C_pilot61, chan[i].prn);
                        codegen_B1C_secondary_code(chan[i].B1C_pilot_sub_code, chan[i].prn);

                        // Generate navigation message
                        // generateNavMsg(grx, &chan[i], advance_fptr);
                        // 根据星历和电离层参数，生成导航消息，存储到 chan[i] 中
                        // generateINavMsg(brx, &chan[i], &eph);
                        generateB1CNavMsg(brx, &chan[i], &eph);

                        // Initialize pseudorange
                        // 计算卫星和接收机之间的伪距
                        computeRange(&rho, eph, brx, xyz, chan[i].prn);
                        // 在这里就已经计算过一一次伪距了
                        chan[i].rho0 = rho;

                        // Initialize carrier phase
                        r_xyz = rho.range;

                        computeRange(&rho, eph, brx, ref, chan[i].prn);
                        r_ref = rho.range;

                        // 得到载波的相位
                        phase_ini = (2.0 * r_ref - r_xyz) / LAMBDA_B1C;
                        chan[i].carr_phase = phase_ini - floor(phase_ini);

                        // // chan[i].code_phase = chan[i].prn * 100 % CA_SEQ_LEN_B1C; // 初始码相位，可以根据PRN号设置一个不同的初始值
                        // chan[i].code_phase = 0.0; // 初始码相位，可以根据PRN号设置一个不同的初始值
                        // chan[i].ibit = chan[i].prn * 10 % NAV_SEQ_LEN_B1C; // 初始比特索引，可以根据PRN号设置一个不同的初始值

                        // chan[i].code_phase = 0.0; // 初始码相位，可以根据PRN号设置一个不同的初始值
                        // chan[i].ibit = chan[i].prn * 10 % NAV_SEQ_LEN_B1C; // 初始比特索引，可以根据PRN号设置一个不同的初始值
                        // chan[i].ibit = 0; // 初始比特索引，可以根据PRN号设置一个不同的初始值



                        // 计算卫星发射 接收机收到的 第一个bit的毫秒时间
                        double ms = (brx.sec - rho.range / SPEED_OF_LIGHT) * 1000.0;
                        if (fmod(ms, 10.0) < 0.0)
                            ms += 10.0;
                        // 令ms对10ms取余，得到码相位，double类型
                        // chan[i].code_phase = (fmod(ms, 10.0) / 10.0) * CA_SEQ_LEN_B1C;
                        // // 计算导航bit的索引，int类型（平均每10ms发送一个symbol）
                        // // chan[i].ibit = int(ceil(fmod(ms, 18000) / 10));
                        // chan[i].ibit = int((fmod(ms, 18000) / 18000) * NAV_SEQ_LEN_B1C);

                        cout << "信号通道分配完成: PRN=" << chan[i].prn << ", 通道号=" << i << endl;
                        // cout << "初始码相位: " << chan[i].code_phase << ", 初始比特索引: " << (fmod(ms, 18000) / 18000) * NAV_SEQ_LEN_B1C << endl;
                        // cout << "码相位比例(ms): " << fmod(ms, 10.0) / 10.0 << endl;
                        // cout << "信号发射时间(ms): " << ms << endl;
                        // cout << "传播时延(s): " << rho.range / SPEED_OF_LIGHT << endl;

                        // computeCodePhase(&chan[i], rho, 0.1, brx);

                        // fprintf(stderr, "--%02d %6.1f %5.1f %11.1f %5.5f\n", chan[i].prn, chan[i].azel[0] * R2D, chan[i].azel[1] * R2D, chan[i].rho0.range, brx.sec);
                        // print_eph2(&eph, sv+1);
                        break;
                    }
                }

                // Set satellite allocation channel
                if (i < MAX_CHAN)
                    allocatedSat[sv] = i;
            }
        }
        else if (allocatedSat[sv] >= 0) // Not visible but allocated
        {
            // Clear channel
            chan[allocatedSat[sv]].prn = 0;

            // Clear satellite allocation flag
            allocatedSat[sv] = -1;
        }
    }
    // advance_fptr = true;
    // cout << "Channel Allocation Completed. Total Visible Satellites: " << nsat << endl;
    return (nsat);
}

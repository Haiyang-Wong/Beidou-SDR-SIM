/**
 * \file bds_sdr.cpp
 * \brief BDS SDR仿真主任务
 * \author LackWood Du
 * \date 2025-12-23
 */

#include "bds_sdr.h"
#pragma message("Compiling bds_sdr.cpp")

void *bds_task(sim_t *sim)
{

    cout << "BDS SDR Simulation Task Started." << endl;

    int iduration;                             // 仿真持续时间，单位为1秒
    double llh[3];                             // 用户初始位置：经度、纬度、高度，单位为度和米，信号生成的目标位置
    double llh_new[3];                         // 用户新位置：经度、纬度、高度，单位为度和米，信号生成的目标位置
    vector<ephem_t> eph_vector[MAX_SAT];       // 存储所有卫星的星历数据
    double xyz[USER_MOTION_SIZE][3];           // 存储用户轨迹数据，每行表示一个时刻的ECEF坐标
    int sv;                                    // 卫星prn编号
    ephem_t eph;                               // 临时存储星历数据
    bdstime_t bdt_min, bdt_max;                // 临时存储BDT时间
    datetime_t t0, tmin, tmax;                 // 临时存储日期时间数据
    datetime_t tl;                             // 临时存储日期时间数据
    bdstime_t g0;                              // 仿真起始时间
    bdstime_t brx;                             // 接收机接收信号的时间
    vector<int> current_eph_index(MAX_SAT, 0); // 当前使用的星历索引
    vector<int> old_eph_index(MAX_SAT, 0);     // 上一时刻使用的星历索引
    // double dt = 0.10000002314200000;           // 时间增量，单位为秒
    double dt = 0.1;           // 时间增量，单位为秒
    channel_t channels[MAX_CHAN];              // 信道数组
    double elvmask = 10;                       // 仰角掩模，单位为度
    vector<int> allocated_channels(MAX_SAT);   // 已分配信道的卫星PRN列表
    float cosPhase[1];                         // 当前载波相位余弦值
    float sinPhase[1];                         // 当前载波相位正弦值
    double temp_phase;                         // 临时载波相位变量
    int iq_buff_size;                          // IQ缓冲区大小
    double delt_between_samples;               // 采样点间的相位增量
    short *iq_buff = NULL;                     // I/Q缓冲区指针
    FILE *fp;                                  // 标准C文件指针，用于文件输入输出操作
    char outfile[MAX_CHAR];
    int num_itime; // 当前仿真时间步数

    iduration = sim->opt.iduration;
    llh[0] = sim->opt.llh[0];
    llh[1] = sim->opt.llh[1];
    llh[2] = sim->opt.llh[2];

    // 读取星历文件
    int eph_count = readBdsB1CEphemerisCpp(eph_vector, sim->opt.navfile);

    // 将经纬度坐标转为对应的ECEF坐标
    llh[0] /= R2D;        // 转为弧度
    llh[1] /= R2D;        // 转为弧度
    llh2xyz(llh, xyz[0]); // 这里是最初的轨迹点

    // 终端打印llh坐标和ECEF坐标
    std::cerr << "xyz = "
              << std::setw(11) << std::fixed << std::setprecision(1) << xyz[0][0] << ", "
              << std::setw(11) << xyz[0][1] << ", "
              << std::setw(11) << xyz[0][2] << '\n';
    std::cerr << "llh = "
              << std::setw(11) << std::setprecision(6) << llh[0] * R2D << ", "
              << std::setw(11) << llh[1] * R2D << ", "
              << std::setw(11) << std::setprecision(1) << llh[2] << '\n';

    // 计算星历的可信时间范围，tmin和tmax表示在这段时间范围内所有的星历都是有效的
    compute_bdt_min_max(eph_vector, &bdt_min, &bdt_max, &tmin, &tmax);

    // 检查并设置仿真起始时间
    set_scenario_start_time_bds(&sim->opt.g0, bdt_min, bdt_max, &t0, &tmin, &tmax, sim->opt.timeoverwrite, eph_count, eph_vector);
    g0 = sim->opt.g0;
    cout << "仿真持续时间：" << (double)iduration << " 秒。" << endl;

    brx = g0; // 接收机接收信号的时间从仿真起始时间开始
    // brx.sec = 2.703420840655015e+05 - 17.79;; // 这里硬编码一个时间，用于测试
    // brx.sec = 359298;

    brx = incBdsTime(brx, dt); // 时间前进0.1秒，准备进入主循环

    // 初始化通道
    init_channel(channels, allocated_channels);

    // cout << "eph_count: " << eph_count << endl;
    // 分配通道
    allocateChannel(channels, eph_vector, brx, xyz[0], elvmask, nullptr, current_eph_index, allocated_channels);

    iq_buff_size = NUM_IQ_SAMPLES;                      // 就是0.1s内生成的采样点的数量
    delt_between_samples = 1.0 / (double)TX_SAMPLERATE; // 采样点间的时间间隔，单位为秒

    // Allocate I/Q buffer
    iq_buff = (short *)calloc(2 * iq_buff_size, sizeof(short));

    // 输出文件路径，即生成的 Galileo 信号 I/Q 数据流文件
    strcpy(outfile, sim->opt.outfile);

    // Open output file
    // "-" can be used as name for stdout
    // 判断究竟是将数据写入到文件还是输出到标准输出(stdout)
    if (strcmp("-", outfile))
    {
        // 使用fopen打开文件，"w"表示以写入模式打开文件，"b"表示以二进制模式打开文件
        if (NULL == (fp = fopen(outfile, "wb")))
        {
            fprintf(stderr, "ERROR: Failed to open output file.\n");
            exit(1);
        }
    }
    else
    {
        fp = stdout;
    }

    // ------------------------------- 主仿真循环 -------------------------------
    // itime以0.1秒为步进，遍历整个仿真时间
    num_itime = iduration * 10;
    brx = incBdsTime(brx, dt);
    vector<int> nav_idx(MAX_CHAN, 0); // 导航电文索引
    for (int itime = 1; itime < num_itime; itime++)
    {
        cout << "仿真时间：" << itime * 0.1 << " 秒。" << endl;
        // 更新接收机的位置，没有实现，这里只是假装一下，实际上目标位置还是静止的
        llh_new[0] = sim->opt.llh[0];
        llh_new[1] = sim->opt.llh[1];
        llh_new[2] = sim->opt.llh[2];

        llh[0] = llh_new[0] / R2D; // 转为弧度
        llh[1] = llh_new[1] / R2D; // 转为弧度
        llh[2] = llh_new[2];       // 高度单位为米
        // 将经纬高坐标转为对应的ECEF坐标
        llh2xyz(llh, xyz[itime]);
        // 遍历所有通道，并找到已分配通道的卫星
        // cout << "1111111111111111" << endl;
        for (int i = 0; i < MAX_CHAN; i++)
        {
            // 如果当前通道已经分配了卫星
            if (channels[i].prn > 0)
            {
                range_t rho;
                sv = channels[i].prn - 1;
                eph = eph_vector[sv][current_eph_index[sv]];
                // 计算伪距
                computeRange(&rho, eph, brx, xyz[itime], channels[i].prn);
                channels[i].azel[0] = rho.azel[0]; // 方位角
                channels[i].azel[1] = rho.azel[1]; // 高度角

                // 计算码相位
                computeCodePhase(&channels[i], rho, dt, brx);   
            }
        }
        // if (itime == 10)
        //     return nullptr;

        // cout << "2222222222222222" << endl;
        // iq_buff_size == 3069000，表示100ms内的采样点数量
        // 接收方也一定是按照这个采样率来采样的
        for (int isample = 0; isample < iq_buff_size; isample++)
        {
            // 当前采样点I/Q数值
            int i_qua = 0;
            int q_qua = 0;
            for (int i_chan = 0; i_chan < MAX_CHAN; i_chan++)
            {
                if (channels[i_chan].prn > 0)
                {
                    // 确保主码码片的索引在合法范围内
                    if (channels[i_chan].code_phase >= CA_SEQ_LEN_B1C)
                    {
                        channels[i_chan].code_phase -= CA_SEQ_LEN_B1C;
                        channels[i_chan].ibit++;                  // 主码的每个码片的持续时间是10ms
                        if (channels[i_chan].ibit >= NAV_SEQ_LEN_B1C) // B1C导航电文长度为1800位
                        {
                            channels[i_chan].ibit = 0;

                            sv = channels[i_chan].prn - 1;
                            eph = eph_vector[sv][current_eph_index[sv]];
                            generateB1CNavMsg(brx, &channels[i_chan], &eph);
                        }
                    }
                    // 使用查表法，计算当前载波相位对应的正弦和余弦值
                    int cosPh = cosTable512[((int)(511 * channels[i_chan].carr_phase)) & 511];
                    int sinPh = sinTable512[((int)(511 * channels[i_chan].carr_phase)) & 511];

                    // 计算当前采样点的主码索引
                    //......
                    int code_index_BOC11 = (int)(channels[i_chan].code_phase * 2);
                    int code_index_BOC61 = (int)(channels[i_chan].code_phase * 12);

                    int B1C_data_chip = channels[i_chan].ca_B1C_data[code_index_BOC11];
                    int B1C_pilot11_chip = channels[i_chan].ca_B1C_pilot11[code_index_BOC11];
                    int B1C_pilot61_chip = channels[i_chan].ca_B1C_pilot61[code_index_BOC61];

                    // 计算当前采样点的导航电文数据位
                    int data_bit = channels[i_chan].nav_bit[channels[i_chan].ibit];

                    // 计算当前采样点的次级码
                    int secondary_code = channels[i_chan].B1C_pilot_sub_code[channels[i_chan].ibit];

                    // 计算实部 I 和虚部 Q
                    double A = 0.5 * B1C_data_chip * data_bit + sqrt(1.0 / 11.0) * B1C_pilot61_chip * secondary_code;
                    double B = sqrt(29.0 / 44.0) * B1C_pilot11_chip * secondary_code;

                    // 计算最终的实部和虚部，结合了载波相位的影响
                    double I = (A * cosPh - B * sinPh); // 实部
                    double Q = (A * sinPh + B * cosPh); // 虚部
                    // // 调制，计算I/Q分量
                    // if (i_chan == 0)
                    // {
                    //     i_qua += I;
                    //     q_qua += Q;
                    // }

                    i_qua += I;
                    q_qua += Q;

                    // 更新码相位
                    channels[i_chan].code_phase += channels[i_chan].f_code * delt_between_samples;
                    // 更新载波相位
                    // 这里的f_carr实际上是载波频偏
                    channels[i_chan].carr_phase += channels[i_chan].f_carr * delt_between_samples;
                    channels[i_chan].carr_phase -= (long)(channels[i_chan].carr_phase); // 保持在0-1范围内

                    // cout << "f_code: " << channels[i_chan].f_code << ", f_carr: " << channels[i_chan].f_carr << endl;
                }
            }

            iq_buff[isample * 2] = (short)i_qua;
            iq_buff[isample * 2 + 1] = (short)q_qua;
            // advance_fptr = true;
        }

        fwrite(iq_buff, sizeof(short), 2 * iq_buff_size, fp);


        // 每隔30秒重新分配一次通道
        int ibrx = (int)(brx.sec * 10 + 0.5);
        if ((int)fmodf(ibrx, 300) == 0)
        {
            allocateChannel(channels, eph_vector, brx, xyz[itime], elvmask, nullptr, current_eph_index, allocated_channels);
        }

        // 时间前进0.1秒
        brx = incBdsTime(brx, dt);

    }

    cout << " BDS SDR Simulation Task Finished. " << endl;
    return nullptr;
}

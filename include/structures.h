/*! \file structures.h
 *  \brief GNSS 结构体定义头文件
 *  \author LackWood Du
 *  \date 2025-12-19
 */
#pragma once

#include <string>
#include "bds_sdr.h"

// 对应 B1C 信号参数
struct B1CParams
{
    int prn;
    int w;
    int p;

    std::string head24_oct; // ICD 表中的“头24个码片（八进制）”
    std::string tail24_oct; // ICD 表中的“末尾24个码片（八进制）”
};

struct datetime_t
{
    int y;      /*!< Calendar year */
    int m;      /*!< Calendar month */
    int d;      /*!< Calendar day */
    int hh;     /*!< Calendar hour */
    int mm;     /*!< Calendar minutes */
    double sec; /*!< Calendar seconds */
};

struct bdstime_t
{
    int week;   // BDT周数
    double sec; // 周内秒数
};

struct ephem_t
{
    int vflg; // 数据有效标志。用于指示当前存储的这组星历参数是否完整、有效，可供计算使用。-

    int PRN;      // 卫星伪随机码编号
    datetime_t t; // -
    bdstime_t time;
    bdstime_t toe; // 轨道参数参考时间 -toe，其中包含Toe，即BDT周内秒，BDT Week，即BDT周数
    int iode;      // 星历版本号 -AODE(IODE)

    bdstime_t toc; // 卫星钟差参考时间 -Toc
    // 描述卫星钟差，即卫星载波时钟与系统时钟之间的差异
    double af0; // 钟偏差 [s] -a0
    double af1; // 钟漂移 [s/s] -a1
    double af2; // 钟加速度 [s/s^2] -a2

    // 下面的参数定义了卫星在地球周围的椭圆轨道
    double deltan; /*!< Delta-N (radians/sec) -Delta_n */
    double cuc;    /*!< Cuc (radians) -Cuc */
    double cus;    /*!< Cus (radians) -Cus */
    double cic;    /*!< Correction to inclination cos (radians) -Cic */
    double cis;    /*!< Correction to inclination sin (radians) -Cis */
    double crc;    /*!< Correction to radius cos (meters) -Crc */
    double crs;    /*!< Correction to radius sin (meters) -Crs */
    double ecc;    // 轨道离心率，描述轨道的椭圆形状，值越接近0表示轨道越接近圆形 -e
    double sqrta;  /*!< sqrt(A) (sqrt(m)) -sqrt(A) */
    double m0;     /*!< Mean anamoly (radians) -M0 */
    double omg0;   /*!< Longitude of the ascending node (radians) -OMEGA0 */
    double inc0;   /*!< Inclination (radians) -i0 */
    double aop;    // 这里代表近地点幅角 -omega
    double omgdot; /*!< Omega dot (radians/s) -OMEGA DOT */
    double idot;   /*!< IDOT (radians/s) -IDOT */

    short ura;  // 描述该卫星星历精度估计值 -SV accuracy
    int svhlth; // 卫星健康状况（0=健康，可用；非0=故障）-SatH1

    // BDS群延迟参数
    double tgd1; // B1/B3频段群延迟 -
    double tgd2; // B2/B3频段群延迟 -

    double gps_time; // -Transmission time，导航电文发送给时刻的BDT时间，该字段由接收机计算并填充得到，非卫星导航电文原始数据
    // 时钟版本号和数据龄
    int IODC; // -AODC

    short flag;

    // Galileo遗留参数
    int codeL2;
    // Working variables follow
    // 中间计算结果（避免重复计算）
    double n;       // 平均运动（平均角速度）
    double sq1e2;   /*!< sqrt(1-e^2) */
    double A;       // 半长轴 A
    double omgkdot; // 相对地球自转速率
    double omg_t;   // 当前时刻的升交点赤经
    int svid;       // 卫星ID

    int satype; // 卫星类型

    // // 不同频率信号传播时延略有不同，这些参数用于电离层与频率间偏差修正
    // double bgde5a; // E5a频段群延迟（Galileo独有）
    // double bgde5b; // E5b频段群延迟
    // double tgd_ext[5]; // 可能扩展的TGD项
};

struct range_t
{
    bdstime_t g;
    double range;      // 伪距
    double rate;       // 伪距变化率 (pseudorange rate)
    double d;          // 卫星与接收机之间的几何距离
    double azel[2];    // 方位角与高度角 (azimuth & elevation)
    double iono_delay; // 电离层延迟
};

struct channel_t
{
    int prn;                   // 伪随机码
    short *ca_B1C_data;        // 指向B1C数据分量伪随机码序列的指针
    short *ca_B1C_pilot11;     // 指向B1C导频分量伪随机码序列的指针
    short *ca_B1C_pilot61;     // 指向B1C导频分量伪随机码序列的指针
    short *B1C_pilot_sub_code; // 指向次级码序列的指针
    double f_carr;             // 载波频率，表示接收机本地生成的用于跟踪卫星信号载波的频率
    double f_code;             // 频码率，表示伪随机码的速率，定义了每个码片的持续时间：T_chip = 1 / f_code
    short *nav_bit;            // 它是指向导航电文页缓冲区的指针
    double carr_phase;         // 载波相位
    double code_phase_p;       // 主码相位
    double code_phase_s;       // 辅助码相位
    bdstime_t g0;              /*!< GPS time at start */
    int ipage;                 // 初始页面索引。导航电文解析的起始页面
    int ibit;                  // 初始比特索引。导航电文解析的起始比特
    int icode;                 // 初始码片索引。导航电文解析的起始码片
    int dataBit;               // 当前数据比特。当前正在处理的导航数据比特值
    float B1C_data_chip;       // 当前E1B码片值。当前时刻E1B信号分量的伪随机码片值。
    float B1C_pilot_chip;      // 当前E1C码片值。当前时刻E1C信号分量的伪随机码片值。
    double azel[2];            // 方位角和仰角。通常 azel[0]为方位角，azel[1]为仰角，表示卫星在天空中的位置。
    range_t rho0;
    bool set_code_phase; // 码相位设置标志。一个布尔标志，可能用于指示是否需要对码相位进行初始化或重置。
    double code_phase;   // ​码相位。伪随机码的相位值。
};

/**
 * @brief 它是GNSS信号仿真器的配置参数容器
 * @note 具体来讲，它告诉程序从哪里读取星历，仿真哪段时间、哪个位置，用什么模式运行，输出到哪里，以及是否通过USRP发射
 */
struct option_t
{
    char navfile[MAX_CHAR]; // 导航星历文件路径
    char umfile[MAX_CHAR];  // 用户模型文件路径
    char outfile[MAX_CHAR]; // 输出文件路径
    char tvfile[MAX_CHAR];  // 时间变量文件路径

    int staticLocationMode; // 是否使用静态位置模式
    int nmeaGGA;            // 是否输出NMEA GGA格式的位置信息
    double iduration;          // 仿真持续时间，单位为1秒
    int verb;               // 是否输出详细日志
    bdstime_t g0;           // 仿真起始时间，定义模拟信号的起始历元
    double llh[3];          // 用户初始位置：经度、纬度、高度，单位为度和米，信号生成的参考位置
    int interactive;        // 是否启用交互模式，允许用户在仿真过程中动态调整参数
    int timeoverwrite;      // 是否符覆盖输入文件中的时间信息，用于强制设仿真时间
    int iono_enable;        // 是否启用电离层延迟模型
    bool use_usrp;          // 是否通过USRP设备发射信号
    bool use_bit_stream;    // 是否使用比特流输入而非浮点信号
};

/**
 * @brief 用于管理 USRP 发射线程的配置和状态，包括线程句柄、互斥锁、UHD 元数据、数据流对象和发送缓冲区。
 * @note tx表示transmit（发射），所以这个结构体主要用于处理与信号发射相关的任务。
 */
struct tx_t
{
    thread worker; // 线程句柄，它是发射线程的标识符
    mutex lock;    // 线程互斥锁，用于保护共享资源，防止多个线程同时访问导致数据不一致
    // int error;
    // uhd::tx_metadata_t md; // UHD（USRP硬件驱动）传输元数据，包含发送数据时的相关信息
    // uhd::tx_streamer::sptr stream; // UHD传输流对象，用于管理数据的发送
    short *buffer;           // 发送缓冲区，存储要发送的数据样本
    const void **buffer_ptr; // 指向发送缓冲区的指针数组，用于传递给UHD发送函数
};

struct bds_t
{
    thread worker;
    mutex lock;
    // int error;

    int ready;
    // pthread_cond_t initialization_done;
};

struct sim_t
{
    option_t opt; // 仿真选项和参数

    tx_t tx;       // USRP发射器配置和状态
    bds_t bds_sim; // BDS信号仿真线程管理

    int status;           // 仿真系统当前状态码，(可能表示运行、暂停、错误等状态)
    short udp_port;       // 用于外部通信的UDP端口号
    bool finished;        // 指示仿真是否已经完成
    short *fifo;          // FIFO环形缓冲区，用于在信号生成线程和USRP发射线程之间传递数据
    long head, tail;      // FIFO缓冲区的读写指针位置
    size_t sample_length; // 表示每一批样本的长度

    condition_variable fifo_read_ready;  // 当FIFO有数据时可通知发射线程
    condition_variable fifo_write_ready; // 当FIFO有空间时可通知写入线程

    double time; // 当前仿真时间，单位为秒
};
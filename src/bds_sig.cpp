
/**
 * \file bds_sig.cpp
 * \brief BDS 伪距和码相位
 * \author LackWood Du
 * \date 2025-12-19
 */

#include "bds_sdr.h"
#pragma message("Compiling bds_sig.cpp")

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
    // pos[3]卫星在某一时刻的位置、速度和钟差ECEF坐标
    // vel[3]卫星速度分量
    // clk[2]卫星是时钟偏差
    double pos[3] = {0.0}, vel[3] = {0.0}, clk[2] = {0.0};
    double los[3] = {0.0};
    double tau = 0.0;
    double range = 0.0, rate = 0.0;
    double xrot = 0.0, yrot = 0.0;

    double llh[3] = {0.0}, neu[3] = {0.0};
    double tmat[3][3] = {0};

    // SV position at time of the pseudorange observation.
    // 获取卫星的位置、速度、钟差
    satpos(eph, g, pos, vel, clk);

    // cout << "PRN: " << prn << ", 卫星位置: "
    //      << std::setw(11) << std::fixed << std::setprecision(1) << pos[0] << ", "
    //      << std::setw(11) << pos[1] << ", "
    //      << std::setw(11) << pos[2] << endl;

    // Receiver to satellite vector and light-time.
    subVect(los, pos, xyz);
    // 向量距离的模长与光速的比值，得到信号传播时间
    tau = normVect(los) / SPEED_OF_LIGHT;

    // Extrapolate the satellite position backwards to the transmission time.
    // pos[0-2]计算的是当前时刻卫星的位置，但是并不是卫星发射信号时的位置，因此需要考虑tau时间内卫星的位移，进行修正。
    // 这里存在问题，按照代码逻辑，tau 是基于卫星在“接收时刻”位置计算的，实际上我们真正需要的是光从卫星发射到接收机的传播时间，也就是信号在卫星发射时刻的位置到接收机的距离除以光速。
    // 因此这种修正方式，还是会存在一定的误差
    pos[0] -= vel[0] * tau;
    pos[1] -= vel[1] * tau;
    pos[2] -= vel[2] * tau;

    // Earth rotation correction. The change in velocity can be neglected.
    // 修正因地球自转带来的误差
    // 这里同样是如此，tau并不是实际误差时间，只是一个近似值，同上
    xrot = pos[0] + pos[1] * GNSS_OMEGA_EARTH_DOT * tau;
    yrot = pos[1] - pos[0] * GNSS_OMEGA_EARTH_DOT * tau;
    pos[0] = xrot;
    pos[1] = yrot;

    // New observer to satellite vector and satellite range.
    // 重新计算接收机到卫星的向量和距离
    subVect(los, pos, xyz);
    range = normVect(los);

    // range = sqrt(std::pow((xyz[0] - pos[0]), 2) + std::pow((xyz[1] - pos[1]),
    // 2) + std::pow((xyz[2] - pos[2]), 2));

    // 卫星和接收机之间的几何距离
    // 虽然上述计算方式存在误差，但是可接受，后面的计算中，就将这个几何距离当作真实距离来使用
    rho->d = range;

    // Pseudorange.
    // 这里是已经知道接收机和卫星之间的实际距离来反推伪距，目的是为了得到伪距之后，生成对应的码片相位信息
    // 这里是将因时钟误差而产生的距离误差加入到几何距离中，来模拟实际卫星时钟误差的影响。
    rho->range = range - SPEED_OF_LIGHT * clk[0];

    double r[3] = {range, range - SPEED_OF_LIGHT * clk[0], 0};

    // Azimuth and elevation angles.
    double satLLH[3];
    xyz2llh(xyz, llh);    // convert userXYZ to llh
    xyz2llh(pos, satLLH); // convert satXYZ to llh
    ltcmat(llh, tmat);
    ecef2neu(los, tmat, neu);
    neu2azel(rho->azel, neu);

    // // Add ionospheric delay
    // // 添加电离层延迟
    // double frequency = CARR_FREQ;
    // rho->iono_delay = ionosphericDelay(ionoutc, g, llh, satLLH, rho->azel, frequency);

    // 为伪距添加电离层延迟参数
    // rho->range += rho->iono_delay;
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
    bdstime_t grx) // checked
{
    double ms;
    int ims;
    double rhorate;

    // chan->rho0.range = rho1.range; // 测试函数功能，后面可以删除这句代码

    // Pseudorange rate.
    // 计算伪距变化率，当卫星靠近接收机时，速度为负值，计算出的频率将增加
    rhorate = (rho1.range - chan->rho0.range) / dt;

    // cout << "PRN: " << chan->prn
    //      << ", 伪距: " << rho1.range
    //      << ", 变化率: " << rhorate
    //      << ", 多普勒频率: " << -rhorate / LAMBDA_B1C << " Hz" << endl;
    // Carrier and code frequency.
    // 载波频偏
    chan->f_carr = (-rhorate / LAMBDA_B1C); // + GALILEO_E1_SUB_CARRIER_A_RATE_HZ;
    // chan->f_carr = 0.0; // 测试时暂时屏蔽多普勒频偏的影响

    // 添加多普勒频偏对码频率的影响，根据载波频偏按比例更改码频率
    // chan->f_code = CODE_FREQ_B1C; // + chan->f_carr * CARR_TO_CODE_B1C;
    chan->f_code = CODE_FREQ_B1C + chan->f_carr * CARR_TO_CODE_B1C;

    // 计算卫星发射 接收机收到的 第一个bit的毫秒时间
    ms = (grx.sec - rho1.range / SPEED_OF_LIGHT) * 1000.0;

    if (fmod(ms, 10.0) < 0.0)
        ms += 10.0;

    // // // 令ms对10ms取余，得到码相位，double类型
    chan->code_phase = (fmod(ms, 10.0) / 10.0) * CA_SEQ_LEN_B1C;

    // 计算导航bit的索引，int类型（平均每10ms发送一个symbol）
    // chan->ibit = int(ceil(fmod(ms, 18000) / 10));
    chan->ibit = int((fmod(ms, 18000) / 18000) * NAV_SEQ_LEN_B1C);
    



    // // Save current pseudorange
    chan->g0 = grx;

    // 修正上一时刻的伪距为当前伪距，为下一次计算做准备
    chan->rho0 = rho1;
    return;
}

/**
 * \brief 将十六进制字符串转换为B1C码片序列
 * \param[out] tmp_ca 存储转换后B1C码片序列
 * \param[in] prn 卫星PRN号
 * \param[in] flag 指定转换的码类型（主码或辅码），取值为 B1C_DATA_PRIMARY、B1C_PILOT_PRIMARY、B1C_PILOT_SUB
 */
void hex_to_b1c_ca(
    std::vector<short> &tmp_ca,
    int prn,
    int flag)
{
    const std::string *hex_str = nullptr;
    int target_len = 0;

    // 1. 选择数据源 + 目标码长
    switch (flag)
    {
    case B1C_DATA_PRIMARY:
        hex_str = &B1C_data_primary_hex_ca.at(prn - 1);
        target_len = CA_SEQ_LEN_B1C;
        break;

    case B1C_PILOT_PRIMARY:
        hex_str = &B1C_pilot_primary_hex_ca.at(prn - 1);
        target_len = CA_SEQ_LEN_B1C;
        break;

    case B1C_PILOT_SUB:
        hex_str = &B1C_pilot_sub_hex_ca.at(prn - 1);
        target_len = SC_SEQ_LEN_B1C;
        break;

    default:
        throw std::invalid_argument("Invalid B1C flag");
    }

    tmp_ca.clear();
    tmp_ca.reserve(target_len);

    // 2. hex → bit → ±1（0 → +1, 1 → -1）
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

        // MSB → LSB
        for (int i = 3; i >= 0; --i)
        {
            if (tmp_ca.size() >= static_cast<size_t>(target_len))
            {
                break;
            }

            int bit = (value >> i) & 0x1;

            // 核心修改点
            tmp_ca.push_back(bit ? -1 : 1);
        }
    }
}

/**
 * \brief 将B1C码片序列转换为BOC(1,1)格式
 * \param[in] tmp_ca 输入的B1C码片序列
 * \param[out] ca 输出的BOC(1,1)格式码片序列
 */
void BOC11(
    const std::vector<short> &tmp_ca,
    short *ca)
{
    int jj = 0;

    for (size_t ii = 0; ii < tmp_ca.size(); ++ii)
    {
        ca[jj] = -tmp_ca[ii];
        ca[jj + 1] = tmp_ca[ii];
        jj += 2;
    }
}

/**
 * \brief 将B1C码片序列转换为BOC(6,1)格式
 * \param[in] tmp_ca 输入的B1C码片序列
 * \param[out] ca 输出的BOC(6,1)格式码片序列
 */
void BOC61(const std::vector<short> &tmp_ca, short *ca)
{
    // 每个主码 chip 展开为 12 个子码片
    const int BOC61_FACTOR = 12;

    for (size_t jj = 0; jj < tmp_ca.size(); ++jj)
    {
        for (int ii = 0; ii < BOC61_FACTOR; ++ii)
        {
            // 符号交替: 奇数位取 -Primary, 偶数位取 +Primary
            ca[jj * BOC61_FACTOR + ii] = ((ii % 2 == 0) ? -1 : 1) * tmp_ca[jj];
        }
    }
}

void BOCmn(
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

void compareB1C_PRN(int prn, const vector<short> &tmp_ca, const string &csvFile)
{
    ifstream infile(csvFile);
    if (!infile.is_open())
    {
        cerr << "无法打开 CSV 文件: " << csvFile << endl;
        return;
    }

    string line;
    getline(infile, line); // 跳过第一行列名

    bool found = false;
    string prn_str;
    string code_str;

    while (getline(infile, line))
    {
        stringstream ss(line);
        string token;

        // 读取第一列 PRN
        getline(ss, token, ',');
        int csv_prn = stoi(token);

        if (csv_prn == prn)
        {
            found = true;
            // 读取第二列 code
            getline(ss, token, ',');
            code_str = token;
            break;
        }
    }

    infile.close();

    if (!found)
    {
        cerr << "CSV 中没有找到 PRN=" << prn << endl;
        return;
    }

    // 转换 CSV 码为 -1 / 1
    vector<short> csv_code;
    for (char c : code_str)
    {
        if (c == '0')
            csv_code.push_back(-1);
        else if (c == '1')
            csv_code.push_back(1);
        else
        {
            cerr << "CSV 码格式错误: " << c << endl;
            return;
        }
    }

    // 长度检查
    if (csv_code.size() != tmp_ca.size())
    {
        cerr << "长度不一致: CSV长度=" << csv_code.size()
             << ", tmp_ca长度=" << tmp_ca.size() << endl;
        return;
    }

    // 比较
    bool match = true;
    for (size_t i = 0; i < tmp_ca.size(); ++i)
    {
        if (tmp_ca[i] != csv_code[i])
        {
            match = false;
            break;
        }
    }

    if (match)
        cout << "PRN " << prn << ": 一致" << endl;
    else
        cout << "PRN " << prn << ": 不一致" << endl;
}

void compareB1C_BOC_code(
    int prn,
    const short *ca,
    int ca_len,
    const string &csv_path)
{
    ifstream file(csv_path);
    if (!file.is_open())
    {
        cerr << "无法打开 CSV 文件: " << csv_path << endl;
        return;
    }

    string line;

    // 跳过表头
    getline(file, line);

    bool found = false;
    vector<short> ref_code;

    while (getline(file, line))
    {
        stringstream ss(line);
        string prn_str, code_str;

        // 读取 PRN
        getline(ss, prn_str, ',');

        int csv_prn = stoi(prn_str);
        if (csv_prn != prn)
            continue;

        // 读取码序列（带引号）
        getline(ss, code_str);

        // 去掉首尾引号
        if (!code_str.empty() && code_str.front() == '"')
            code_str.erase(0, 1);
        if (!code_str.empty() && code_str.back() == '"')
            code_str.pop_back();

        // 解析 ±1
        stringstream code_ss(code_str);
        int val;
        while (code_ss >> val)
        {
            ref_code.push_back(static_cast<short>(val));
        }

        found = true;
        break;
    }

    file.close();

    if (!found)
    {
        cout << "PRN " << prn << " : CSV 中未找到该 PRN" << endl;
        return;
    }

    if ((int)ref_code.size() != ca_len)
    {
        cout << "PRN " << prn << " : 不一致（长度不匹配，CSV="
             << ref_code.size() << ", 本地=" << ca_len << "）" << endl;
        return;
    }

    for (int i = 0; i < ca_len; ++i)
    {
        if (ref_code[i] != ca[i])
        {
            cout << "PRN " << prn << " : 不一致（索引 " << i
                 << ", CSV=" << ref_code[i]
                 << ", 本地=" << ca[i] << "）" << endl;
            return;
        }
    }

    cout << "PRN " << prn << " : 一致" << endl;
}

void codegen_B1C_data(short *ca, int prn)
{
    vector<short> tmp_ca;
    hex_to_b1c_ca(tmp_ca, prn, B1C_DATA_PRIMARY);
    // compareB1C_PRN(prn, tmp_ca, "D:\\beidou\\BDS_SDR_SIM\\test_data\\Primary_code_parameter\\B1C_Data_Primary_Code.csv");
    // BOC11(tmp_ca, ca);
    BOCmn(tmp_ca, ca, 1, 1);
    // BOCmn(tmp_ca, ca, 6, 1);
    // compareB1C_BOC_code(prn, ca, CA_SEQ_LEN_B1C * 2, "D:\\beidou\\BDS_SDR_SIM\\data\\B1C_BOC11_Data.csv");
}

void codegen_B1C_pilot(short *ca11, short *ca61, int prn)
{
    vector<short> tmp_ca;
    hex_to_b1c_ca(tmp_ca, prn, B1C_PILOT_PRIMARY);

    // compareB1C_PRN(prn, tmp_ca, "D:\\beidou\\BDS_SDR_SIM\\test_data\\Primary_code_parameter\\B1C_Pilot_Primary_Code.csv");
    // BOC11(tmp_ca, ca11);
    BOCmn(tmp_ca, ca11, 1, 1);
    // BOCmn(tmp_ca, ca11, 6, 1);
    // BOC61(tmp_ca, ca61);
    BOCmn(tmp_ca, ca61, 6, 1);
    // BOCmn(tmp_ca, ca61, 1, 1);
}

void codegen_B1C_secondary_code(short *sc, int prn)
{
    vector<short> tmp_ca;
    hex_to_b1c_ca(tmp_ca, prn, B1C_PILOT_SUB);
    // compareB1C_PRN(prn, tmp_ca, "D:\\beidou\\BDS_SDR_SIM\\test_data\\Primary_code_parameter\\B1C_Pilot_Sub_Code.csv");
    // 直接复制次级码
    for (size_t i = 0; i < tmp_ca.size(); ++i)
    {
        sc[i] = tmp_ca[i];
    }
}
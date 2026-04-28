/**
 * @file channel.cpp
 * @brief GNSS 信道处理相关函数实现
 * @author LackWood Du
 * @date 2025-12-24
 */
#include "bds_sdr.h"

int generateINavMsg(bdstime_t g, channel_t *chan, ephem_t *eph, int idx)
{
    cout<< "生成 INav 电文: PRN=" << (eph ? eph->PRN : -1) << ", idx=" << idx << endl;
    if (chan == nullptr || eph == nullptr)
    {
        return -1;
    }

    const int group_count = static_cast<int>(nav_data_18_414M.size());
    if (group_count <= 0)
    {
        return -1;
    }

    // idx 超出范围时按循环方式回绕到 [0, group_count)
    const int wrapped_idx = ((idx % group_count) + group_count) % group_count;

    // 1. 找到对应 PRN 的导航电文
    string nav_bits_str;
    bool found = false;
    for (const auto &p : nav_data_18_414M[wrapped_idx])
    {
        if (p.first == eph->PRN)
        {
            nav_bits_str = p.second;
            found = true;
            break;
        }
    }

    if (!found)
    {
        return -1; // 没有找到对应 PRN
    }

    // 2. 确保 chan->nav_bit 已经分配足够空间
    //    这里假设已经在外部分配了 nav_bit
    //    如果没有，可以在这里分配: chan->nav_bit = new int[nav_bits_str.size()];

    // 3. 将字符 '0'/'1' 转换为 int 存入 chan->nav_bit
    for (size_t i = 0; i < nav_bits_str.size(); i++)
    {
        if (nav_bits_str[i] == '1')
            chan->nav_bit[i] = -1;
        else
            chan->nav_bit[i] = 1; // 如果需要 -1 表示 0，可以改成 -1
    }

    cout << "INav电文生成: PRN=" << eph->PRN << ", idx=" << idx
         << " (wrapped_idx=" << wrapped_idx << ")" << endl;

    return 1; // 成功
}

static void appendBits(short *bits, int &pos, uint64_t value, int width)
{
    for (int i = width - 1; i >= 0; --i)
    {
        bits[pos++] = (short)((value >> i) & 1ULL);
    }
}

static int64_t quantizeSigned(double value, int scale_exp)
{
    // raw = value / 2^scale_exp
    return (int64_t)llround(ldexp(value, -scale_exp));
}

// 如果scale_exp值为负数，表示value将被放大后量化；如果scale_exp值为正数，表示value将被缩小后量化。
// 例如，quantizeSigned(0.5, -3) 将计算 0.5 * 2^3 = 4.0，并返回 4。
static uint64_t quantizeUnsigned(double value, int scale_exp)
{
    double raw = ldexp(value, -scale_exp);
    if (raw < 0.0)
        raw = 0.0;
    return (uint64_t)llround(raw);
}

static uint64_t twosComplement(int64_t value, int width)
{
    if (width >= 64)
        return (uint64_t)value;

    uint64_t mask = (1ULL << width) - 1ULL;
    return ((uint64_t)value) & mask;
}

static uint32_t crc24qBits(const short *bits, int bit_count)
{
    const uint32_t poly = 0x1864CFB;
    uint32_t crc = 0;

    for (int i = 0; i < bit_count; ++i)
    {
        uint32_t bit = (uint32_t)(bits[i] & 1);
        uint32_t top = (crc >> 23) & 1U;

        crc = ((crc << 1) & 0xFFFFFFU);
        if ((top ^ bit) != 0)
            crc ^= poly;
    }

    return crc & 0xFFFFFFU;
}

static int getBdsSatType(const ephem_t *eph)
{
    if (eph->svid <= 5 || eph->svid >= 59)
        return 1; // GEO

    return (eph->A > 4.0e7) ? 2 : 3; // IGSO : MEO
}

static int GF64Mul(int a, int b)
{
    int p = 0;

    a &= 0x3f;
    b &= 0x3f;

    for (int i = 0; i < 6; ++i)
    {
        if (b & 1)
            p ^= a;

        b >>= 1;
        a <<= 1;

        // GF(64), primitive polynomial: x^6 + x + 1
        if (a & 0x40)
            a ^= 0x43;
    }

    return p & 0x3f;
}

// info[0] 是第 1 bit，按 MSB-first 每 6 bit 组成一个 GF(64) 符号。
// coded[0..1199] 输出 200 个 6-bit 符号展开后的 0/1 bit。
void B1CSubframe2LDPCEncode600(const short info[600], short coded[1200])
{
    int symbol[B1C_SF2_CODE_SYMBOLS];
    const char *m = B1CMatrixGen2;

    memset(symbol, 0, sizeof(symbol));

    // 600 bit -> 100 GF(64) information symbols
    for (int i = 0; i < B1C_SF2_INFO_SYMBOLS; ++i)
    {
        int v = 0;
        for (int j = 0; j < 6; ++j)
            v = (v << 1) | (info[i * 6 + j] ? 1 : 0);

        symbol[i] = v;
    }

    // Generate 100 GF(64) parity symbols:
    // parity[i] = sum_j Matrix[i][j] * info_symbol[j] over GF(64)
    for (int i = 0; i < B1C_SF2_INFO_SYMBOLS; ++i)
    {
        int sum = 0;

        for (int j = 0; j < B1C_SF2_INFO_SYMBOLS; ++j)
        {
            int coef = ((int)(*m++)) - '0';
            sum ^= GF64Mul(coef, symbol[j]);
        }

        symbol[B1C_SF2_INFO_SYMBOLS + i] = sum;
    }

    // 200 GF(64) symbols -> 1200 output bits, MSB-first
    for (int i = 0; i < B1C_SF2_CODE_SYMBOLS; ++i)
    {
        for (int j = 0; j < 6; ++j)
            coded[i * 6 + j] = (short)((symbol[i] >> (5 - j)) & 1);
    }
}

static void B1CSubframe3LDPCEncode264(const short info[264], short coded[528])
{
    int symbol[88];
    memset(symbol, 0, sizeof(symbol));

    for (int i = 0; i < 44; ++i)
    {
        int v = 0;
        for (int j = 0; j < 6; ++j)
            v = (v << 1) | (info[i * 6 + j] ? 1 : 0);

        symbol[i] = v;
    }

    // 仿照当前项目 src/BCNav1Bit.cpp 的现有实现：
    // LDPCEncode(Symbol3, B1C_SUBFRAME3_SYMBOL_LENGTH, B1CMatrixGen2);
    const char *m = B1CMatrixGen2;

    for (int i = 0; i < 44; ++i)
    {
        int sum = 0;

        for (int j = 0; j < 44; ++j)
        {
            const int coef = ((int)(*m++)) - '0';
            sum ^= GF64Mul(coef, symbol[j]);
        }

        symbol[44 + i] = sum;
    }

    for (int i = 0; i < 88; ++i)
    {
        for (int j = 0; j < 6; ++j)
            coded[i * 6 + j] = (short)((symbol[i] >> (5 - j)) & 1);
    }
}


int generateB1CSubframe1(bdstime_t g, channel_t *chan, ephem_t *eph, short *subframe)
{
    (void)chan;

    if (subframe == nullptr || eph == nullptr)
    {
        return -1;
    }

    // B1C ICD 规定 subframe 1 里的 PRN 字段长度是 6 bit，有效范围是 1~63， 0x3F == 0011 1111
    const unsigned int prn = static_cast<unsigned int>(eph->PRN) & 0x3Fu;
    const int sec_of_week = static_cast<int>(floor(g.sec));
    int sec_of_hour = sec_of_week % 3600;
    if (sec_of_hour < 0)
    {
        sec_of_hour += 3600;
    }
    const unsigned int soh = static_cast<unsigned int>(sec_of_hour / 18) & 0xFFu;

    // 写指针
    short reg21[6];
    int idx = 5;
    for (int i = 5; i >= 0; --i)
    {
        reg21[idx--] = ((prn >> i) & 1u) ? (short)-1 : (short)1;
    }
    for (int i = 0; i < 21; ++i)
    {
        subframe[i] = static_cast<short>((1 - reg21[5]) >> 1);
        const short feedback = static_cast<short>(reg21[1] * reg21[3] * reg21[4] * reg21[5]);
        reg21[5] = reg21[4];
        reg21[4] = reg21[3];
        reg21[3] = reg21[2];
        reg21[2] = reg21[1];
        reg21[1] = reg21[0];
        reg21[0] = feedback;
    }

    short reg51[8];
    idx = 7;
    for (int i = 7; i >= 0; --i)
    {
        reg51[idx--] = ((soh >> i) & 1u) ? (short)-1 : (short)1;
    }
    for (int i = 0; i < 51; ++i)
    {
        subframe[21 + i] = static_cast<short>((1 - reg51[7]) >> 1);
        const short feedback = static_cast<short>(reg51[0] * reg51[3] * reg51[4] *
                                                 reg51[5] * reg51[6] * reg51[7]);
        reg51[7] = reg51[6];
        reg51[6] = reg51[5];
        reg51[5] = reg51[4];
        reg51[4] = reg51[3];
        reg51[3] = reg51[2];
        reg51[2] = reg51[1];
        reg51[1] = reg51[0];
        reg51[0] = feedback;
    }

    return 1;
}

int generateB1CSubframe2(bdstime_t g, channel_t *chan, ephem_t *eph, short *subframe)
{
    const int subframe_len = 1200;
    const int subframe2_info_len = 600;
    const int subframe2_crc_input_len = 576;

    (void)chan;

    if (subframe == nullptr || eph == nullptr || eph->vflg == 0)
        return -1;

    for (int i = 0; i < subframe_len; i++)
        subframe[i] = 0;

    short info[600];
    memset(info, 0, sizeof(info));

    int pos = 0;

    const int wn = g.week & 0x1FFF;
    const int how = ((int)floor(g.sec / 3600.0)) & 0xFF;
    const int iodc = eph->IODC & 0x3FF;
    const int iode = eph->iode & 0xFF;
    const int toe = ((int)llround(eph->toe.sec / 300.0)) & 0x7FF;
    const int sat_type = getBdsSatType(eph) & 0x3;

    const double axis = (eph->A > 0.0) ? eph->A : eph->sqrta * eph->sqrta;
    const double ref_axis = (sat_type == 3) ? 27906100.0 : 42162200.0;

    appendBits(info, pos, (uint64_t)wn, 13);
    appendBits(info, pos, (uint64_t)how, 8);
    appendBits(info, pos, (uint64_t)iodc, 10);

    appendBits(info, pos, (uint64_t)iode, 8);
    appendBits(info, pos, (uint64_t)toe, 11);
    appendBits(info, pos, (uint64_t)sat_type, 2);

    appendBits(info, pos, twosComplement(quantizeSigned(axis - ref_axis, -9), 26), 26);
    appendBits(info, pos, twosComplement(0, 25), 25); // A_dot

    appendBits(info, pos, twosComplement(quantizeSigned(eph->deltan / PI, -44), 17), 17);
    appendBits(info, pos, twosComplement(0, 23), 23); // Delta_n_dot

    appendBits(info, pos, twosComplement(quantizeSigned(eph->m0 / PI, -32), 33), 33);
    appendBits(info, pos, quantizeUnsigned(eph->ecc, -34), 33);
    appendBits(info, pos, twosComplement(quantizeSigned(eph->aop / PI, -32), 33), 33);

    appendBits(info, pos, twosComplement(quantizeSigned(eph->omg0 / PI, -32), 33), 33);
    appendBits(info, pos, twosComplement(quantizeSigned(eph->inc0 / PI, -32), 33), 33);
    appendBits(info, pos, twosComplement(quantizeSigned(eph->omgdot / PI, -44), 19), 19);
    appendBits(info, pos, twosComplement(quantizeSigned(eph->idot / PI, -44), 15), 15);

    appendBits(info, pos, twosComplement(quantizeSigned(eph->cis, -30), 16), 16);
    appendBits(info, pos, twosComplement(quantizeSigned(eph->cic, -30), 16), 16);
    appendBits(info, pos, twosComplement(quantizeSigned(eph->crs, -8), 24), 24);
    appendBits(info, pos, twosComplement(quantizeSigned(eph->crc, -8), 24), 24);
    appendBits(info, pos, twosComplement(quantizeSigned(eph->cus, -30), 21), 21);
    appendBits(info, pos, twosComplement(quantizeSigned(eph->cuc, -30), 21), 21);

    appendBits(info, pos, (uint64_t)((int)llround(eph->toc.sec / 300.0) & 0x7FF), 11);
    appendBits(info, pos, twosComplement(quantizeSigned(eph->af0, -34), 25), 25);
    appendBits(info, pos, twosComplement(quantizeSigned(eph->af1, -50), 22), 22);
    appendBits(info, pos, twosComplement(quantizeSigned(eph->af2, -66), 11), 11);

    appendBits(info, pos, twosComplement(quantizeSigned(eph->tgd1, -34), 12), 12); // TGD_B2ap, current-project compatible
    appendBits(info, pos, 0, 12); // ISC_B1C                            // ISC_B1C
    appendBits(info, pos, twosComplement(quantizeSigned(eph->tgd1, -34), 12), 12); // TGD_B1C

    appendBits(info, pos, 0, 7); // Rev

    if (pos != 576)
        return -2;

    uint32_t crc = crc24qBits(info, subframe2_crc_input_len);
    appendBits(info, pos, crc, 24);

    if (pos != subframe2_info_len)
        return -3;

    B1CSubframe2LDPCEncode600(info, subframe);
    return 1;

}

int generateB1CSubframe3(bdstime_t g, channel_t *chan, ephem_t *eph, short *subframe)
{
    const int subframe_len = 528;
    const int info_len = 264;
    const int crc_input_len = 240;

    (void)chan;

    if (subframe == nullptr)
        return -1;

    for (int i = 0; i < subframe_len; ++i)
        subframe[i] = 0;

    short info[info_len];
    memset(info, 0, sizeof(info));

    const int frame_index = (int)floor(g.sec / 18.0);
    int soh = frame_index % 200;
    if (soh < 0)
        soh += 200;

    const int index = (soh >= 100) ? (soh - 100) : soh;

    // cout << "生成 B1C Subframe 3: PRN=" << (eph ? eph->PRN : -1) << ", frame_index=" << frame_index
    //      << ", soh=" << soh << ", index=" << index << endl;

    unsigned int page_id = 4;
    if ((index % 10) == 0)
        page_id = 1;
    else if ((index % 10) == 5)
        page_id = 3;
    else if (((index % 5) == 1) && index < 80)
        page_id = 2;
    else if (index == 99)
        page_id = 0;

    unsigned int flags = 0;
    if (eph != nullptr && eph->svhlth != 0)
        flags |= 0x80u; // HS

    int pos = 0;

    // 对齐当前项目 ComposeSubframe3():
    // Frame3Data[0] = Flags << 9;
    // Frame3Data[0] |= PageID << 18;
    appendBits(info, pos, 0, 3);        // word0 bit23..21
    appendBits(info, pos, page_id, 3);  // word0 bit20..18: PageID
    appendBits(info, pos, 0, 1);        // word0 bit17
    appendBits(info, pos, flags, 8);    // word0 bit16..9
    appendBits(info, pos, 0, 9);        // word0 bit8..0

    // 所有类型页面的数据字段全部为 0，补齐前 10 个 24-bit word。
    appendBits(info, pos, 0, crc_input_len - pos);

    if (pos != crc_input_len)
        return -2;

    const uint32_t crc = crc24qBits(info, crc_input_len);
    appendBits(info, pos, crc, 24);

    if (pos != info_len)
        return -3;

    B1CSubframe3LDPCEncode264(info, subframe);
    return 1;
}

static void interleaveB1CNavMsg(
    const short subframe1[72],
    const short subframe2[1200],
    const short subframe3[528],
    short nav_msg[1800])
{
    memset(nav_msg, 0, 1800 * sizeof(short));

    // Subframe 1: first 72 bits, no interleaving.
    for (int i = 0; i < 72; ++i)
        nav_msg[i] = subframe1[i];

    // 对齐当前项目 BCNav1Bit.cpp:
    // 11 round of subframe2, subframe2, subframe3
    for (int i = 0; i < 11; ++i)
    {
        short *p1 = nav_msg + 72 + i * 3;
        short *p2 = nav_msg + 73 + i * 3;
        short *p3 = nav_msg + 74 + i * 3;

        for (int j = 0; j < 48; ++j)
        {
            *p1 = subframe2[i * 96 + j];
            *p2 = subframe2[i * 96 + 48 + j];
            *p3 = subframe3[i * 48 + j];

            p1 += 36;
            p2 += 36;
            p3 += 36;
        }
    }

    // Last three rows of subframe2.
    short *p1 = nav_msg + 105;
    short *p2 = nav_msg + 106;
    short *p3 = nav_msg + 107;

    for (int j = 0; j < 48; ++j)
    {
        *p1 = subframe2[22 * 48 + j];
        *p2 = subframe2[23 * 48 + j];
        *p3 = subframe2[24 * 48 + j];

        p1 += 36;
        p2 += 36;
        p3 += 36;
    }
}


int generateB1CNavMsg(bdstime_t g, channel_t *chan, ephem_t *eph)
{
    const int nav_msg_len = 1800;

    if (chan == nullptr || eph == nullptr || chan->nav_bit == nullptr)
        return -1;

    short subframe1[72];
    short subframe2[1200];
    short subframe3[528];
    short nav_msg[1800];

    if (generateB1CSubframe1(g, chan, eph, subframe1) < 0 ||
        generateB1CSubframe2(g, chan, eph, subframe2) < 0 ||
        generateB1CSubframe3(g, chan, eph, subframe3) < 0)
    {
        return -1;
    }

    interleaveB1CNavMsg(subframe1, subframe2, subframe3, nav_msg);

    // 0/1 bit -> signal symbol: 0 => +1, 1 => -1
    for (int i = 0; i < nav_msg_len; ++i)
        chan->nav_bit[i] = nav_msg[i] ? (short)-1 : (short)1;

    return 1;
}


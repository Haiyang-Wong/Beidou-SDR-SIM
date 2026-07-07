#include "bds_sdr.h"

static constexpr int B1C_SUBFRAME1_BITS = 72;
static constexpr int B1C_SUBFRAME2_DATA_BITS = 576;
static constexpr int B1C_SUBFRAME2_INFO_BITS_LOCAL = 600;
static constexpr int B1C_SUBFRAME2_CODE_BITS_LOCAL = 1200;
static constexpr int B1C_SUBFRAME3_DATA_BITS = 240;
static constexpr int B1C_SUBFRAME3_INFO_BITS = 264;
static constexpr int B1C_SUBFRAME3_CODE_BITS = 528;
static constexpr int B1C_NAV_BITS = 1800;
static constexpr int B1C_ALMANAC_COUNT = 63;
static constexpr double B1C_GEO_IGSO_REF_AXIS = 42162200.0;
static constexpr double B1C_MEO_REF_AXIS = 27906100.0;

static void assembleFrameBits(short *bits, int &pos, uint64_t value, int width)
{
    if (width <= 0)
        return;

    if (width < 64)
        value &= ((1ULL << width) - 1ULL);

    for (int i = width - 1; i >= 0; --i)
    {
        if (i >= 64)
            bits[pos++] = 0;
        else
            bits[pos++] = static_cast<short>((value >> i) & 1ULL);
    }
}

static void assembleSignedFrameBits(short *bits, int &pos, int64_t value, int width)
{
    const uint64_t mask = (width >= 64) ? ~0ULL : ((1ULL << width) - 1ULL);
    assembleFrameBits(bits, pos, static_cast<uint64_t>(value) & mask, width);
}

static void assembleFrameBitBlock(short *dst, int &pos, const short *src, int count)
{
    for (int i = 0; i < count; ++i)
        dst[pos++] = src[i] ? 1 : 0;
}

static uint64_t scaleValueToUnsignedField64(double value, int scale)
{
    const double abs_value = fabs(value);
    if (abs_value == 0.0)
        return 0;

    uint64_t raw = 0;
    memcpy(&raw, &abs_value, sizeof(raw));

    const int exp = static_cast<int>((raw >> 52) & 0x7ffULL);
    uint64_t fraction = raw & 0x000fffffffffffffULL;
    if (exp == 0 && fraction == 0)
        return 0;

    fraction |= 0x0010000000000000ULL;
    const int shift = 1074 - exp + scale;

    if (shift >= 0 && shift < 63)
    {
        fraction += (1ULL << shift);
        fraction >>= (shift + 1);
        return fraction;
    }

    const long double scaled = ldexp(static_cast<long double>(abs_value), -scale);
    if (scaled <= 0.0L)
        return 0;

    return static_cast<uint64_t>(floor(scaled + 0.5L));
}

static int64_t scaleValueToSignedField64(double value, int scale)
{
    const uint64_t fraction = scaleValueToUnsignedField64(value, scale);
    return std::signbit(value) ? -static_cast<int64_t>(fraction) : static_cast<int64_t>(fraction);
}

static int32_t scaleValueToSignedField32(double value, int scale)
{
    return static_cast<int32_t>(scaleValueToSignedField64(value, scale));
}

static uint32_t scaleValueToUnsignedField32(double value, int scale)
{
    return static_cast<uint32_t>(scaleValueToUnsignedField64(value, scale));
}

static uint32_t crc24qBits(const short *bits, int bit_count)
{
    uint32_t crc = 0;
    const uint32_t poly = 0x1864CFBu;

    for (int i = 0; i < bit_count; ++i)
    {
        crc ^= static_cast<uint32_t>(bits[i] & 1) << 23;
        if ((crc & 0x800000u) != 0)
            crc = ((crc << 1) ^ poly) & 0xffffffu;
        else
            crc = (crc << 1) & 0xffffffu;
    }

    return crc & 0xffffffu;
}

static int gf64Mul(int a, int b)
{
    int product = 0;

    a &= 0x3f;
    b &= 0x3f;

    for (int i = 0; i < 6; ++i)
    {
        if ((b & (1 << i)) != 0)
            product ^= a << i;
    }

    for (int k = 10; k >= 6; --k)
    {
        if ((product & (1 << k)) != 0)
        {
            product ^= 1 << k;
            product ^= 1 << (k - 5);
            product ^= 1 << (k - 6);
        }
    }

    return product & 0x3f;
}

static void bitsToSymbols6(const short *bits, int symbol_count, int *symbols)
{
    for (int i = 0; i < symbol_count; ++i)
    {
        int value = 0;
        for (int j = 0; j < 6; ++j)
            value = (value << 1) | (bits[i * 6 + j] ? 1 : 0);
        symbols[i] = value;
    }
}

static void symbols6ToBits(const int *symbols, int symbol_count, short *bits)
{
    int pos = 0;
    for (int i = 0; i < symbol_count; ++i)
    {
        for (int j = 5; j >= 0; --j)
            bits[pos++] = static_cast<short>((symbols[i] >> j) & 1);
    }
}

static void ldpcEncodeSymbols(int *symbols, int info_symbol_count, const char *matrix)
{
    int *parity = symbols + info_symbol_count;

    for (int row = 0; row < info_symbol_count; ++row)
    {
        int sum = 0;
        for (int col = 0; col < info_symbol_count; ++col)
        {
            const int coef = static_cast<int>(matrix[row * info_symbol_count + col]) - '0';
            sum ^= gf64Mul(coef, symbols[col]);
        }
        parity[row] = sum;
    }
}

static uint32_t bchPrn21(int svid)
{
    static const uint32_t basis[6] = {
        0x00a4cbu,
        0x014996u,
        0x0237e7u,
        0x046fceu,
        0x087b57u,
        0x105265u};

    uint32_t value = 0;
    for (int i = 0; i < 6; ++i)
    {
        if ((svid & (1 << i)) != 0)
            value ^= basis[i];
    }

    return value & 0x1fffffu;
}

static uint64_t bchSoh51(int soh)
{
    static const uint64_t basis[8] = {
        0x00f3a905b4be3ULL,
        0x0114fb0eddc25ULL,
        0x0229f61dbb84aULL,
        0x0453ec3b77094ULL,
        0x085471735aacbULL,
        0x105b4be301e75ULL,
        0x20453ec3b7709ULL,
        0x4079d482da5f1ULL};

    uint64_t value = 0;
    for (int i = 0; i < 8; ++i)
    {
        if ((soh & (1 << i)) != 0)
            value ^= basis[i];
    }

    return value & 0x7ffffffffffffULL;
}

static unsigned char getB1CSatType(const ephem_t *eph)
{
    if (eph == nullptr)
        return 0;

    if (eph->flag != 0)
        return static_cast<unsigned char>(eph->flag & 0x3u);

    const double axis = (eph->axis > 0.0) ? eph->axis : eph->sqrt_a * eph->sqrt_a;
    if (eph->svid <= 5 || eph->svid >= 59)
        return 1;
    return (axis > 4.0e7) ? 2 : 3;
}

static void prepareB1CEphemeris(ephem_t &eph)
{
    if (eph.axis <= 0.0 && eph.sqrt_a > 0.0)
        eph.axis = eph.sqrt_a * eph.sqrt_a;

    if (eph.flag == 0)
        eph.flag = getB1CSatType(&eph);

    if (eph.n == 0.0 && eph.sqrt_a > 0.0 && eph.axis > 0.0)
        eph.n = CGCS2000_SQRT_GM / (eph.sqrt_a * eph.axis) + eph.delta_n;

    if (eph.tgd_ext[0] == 0.0 && eph.tgd_ext[1] == 0.0 && eph.tgd != 0.0)
    {
        eph.tgd_ext[0] = eph.tgd;
        eph.tgd_ext[1] = eph.tgd;
    }

    if ((eph.toe % 300) != 0 && eph.sqrt_a > 0.0 && eph.axis > 0.0)
    {
        int new_toe = (eph.toe + 150) / 300 * 300;
        const int time_diff = new_toe - eph.toe;
        double mean_motion = CGCS2000_SQRT_GM / (eph.sqrt_a * eph.axis) + eph.delta_n;

        if (new_toe >= static_cast<int>(SECONDS_IN_WEEK))
        {
            new_toe -= static_cast<int>(SECONDS_IN_WEEK);
            ++eph.week;
        }

        eph.toe = new_toe;
        eph.toc = new_toe;

        eph.axis += eph.axis_dot * time_diff;
        if (eph.axis > 0.0)
            eph.sqrt_a = sqrt(eph.axis);

        eph.delta_n += eph.delta_n_dot * time_diff;
        mean_motion += eph.delta_n_dot * time_diff;
        eph.n = mean_motion;
        eph.m0 += mean_motion * time_diff;
        eph.i0 += eph.idot * time_diff;
        eph.omega0 += eph.omega_dot * time_diff;
    }
}

static uint32_t b1cFlagsFromEph(const ephem_t *eph)
{
    if (eph == nullptr)
        return 0;

    const uint32_t health_flag = (eph->health & 0x20u) ? 0x80u : 0u;
    const uint32_t integrity = eph->flag;
    return health_flag | ((integrity & 0x1cu) << 2) | (integrity >> 11);
}

int encodeB1CSubframe1Bits(const short nav_bits[14], short encoded_bits[72])
{
    if (nav_bits == nullptr || encoded_bits == nullptr)
        return -1;

    int svid = 0;
    int soh = 0;

    for (int i = 0; i < 14; ++i)
    {
        if (nav_bits[i] != 0 && nav_bits[i] != 1)
            return -1;

        if (i < 6)
            svid = (svid << 1) | nav_bits[i];
        else
            soh = (soh << 1) | nav_bits[i];
    }

    int pos = 0;
    assembleFrameBits(encoded_bits, pos, bchPrn21(svid), 21);
    assembleFrameBits(encoded_bits, pos, bchSoh51(soh), 51);
    return (pos == B1C_SUBFRAME1_BITS) ? 1 : -2;
}

b1c_nav_time_fields_t computeB1CNavTimeFields(bdstime_t g)
{
    constexpr double frame_period_ms = 18.0 * BDS_MILLISECONDS_IN_SECOND;
    constexpr double eps_ms = 1e-3;

    normalizeBdsTime(g);

    b1c_nav_time_fields_t fields{};
    fields.tx_ms_raw = static_cast<double>(g.milliseconds) + g.sub_milliseconds;
    fields.tx_sow_raw = fields.tx_ms_raw / BDS_MILLISECONDS_IN_SECOND;

    double normalized_ms = fields.tx_ms_raw + eps_ms;
    if (normalized_ms >= BDS_MILLISECONDS_IN_WEEK)
        normalized_ms -= BDS_MILLISECONDS_IN_WEEK;
    else if (normalized_ms < 0.0)
        normalized_ms += BDS_MILLISECONDS_IN_WEEK;

    fields.tx_ms_normalized = normalized_ms;
    fields.tx_sow_normalized = normalized_ms / BDS_MILLISECONDS_IN_SECOND;

    const int frame_index = static_cast<int>(normalized_ms / frame_period_ms);
    const int how = frame_index / 200;
    const int soh = frame_index % 200;

    fields.how = how & 0xff;
    fields.soh = soh & 0xff;
    fields.floor_result = soh;
    fields.hour_start_ms = static_cast<double>(how) * BDS_MILLISECONDS_IN_HOUR;
    fields.hour_start_sow = fields.hour_start_ms / BDS_MILLISECONDS_IN_SECOND;
    fields.frame_ratio = (normalized_ms - fields.hour_start_ms) / frame_period_ms;
    fields.frame_start_ms = static_cast<double>(frame_index) * frame_period_ms;
    fields.frame_start_sow = fields.frame_start_ms / BDS_MILLISECONDS_IN_SECOND;
    fields.tow_from_how_soh = how * static_cast<int>(SECONDS_IN_HOUR) + soh * 18;

    return fields;
}

static int generateSubframe1Bits(int svid, int soh, short *subframe)
{
    if (subframe == nullptr)
        return -1;

    int pos = 0;
    assembleFrameBits(subframe, pos, bchPrn21(svid & 0x3f), 21);
    assembleFrameBits(subframe, pos, bchSoh51(soh & 0xff), 51);
    return (pos == B1C_SUBFRAME1_BITS) ? 1 : -2;
}

int generateB1CSubframe1(bdstime_t g, channel_t *chan, ephem_t *eph, short *subframe)
{
    (void)chan;

    if (subframe == nullptr || eph == nullptr)
        return -1;

    const b1c_nav_time_fields_t time_fields = computeB1CNavTimeFields(g);
    return generateSubframe1Bits(eph->svid, time_fields.soh, subframe);
}

static void generateEphemeris1Bits(const ephem_t &eph, short out[211])
{
    int pos = 0;
    const double ref_axis = (getB1CSatType(&eph) == 3) ? B1C_MEO_REF_AXIS : B1C_GEO_IGSO_REF_AXIS;

    assembleFrameBits(out, pos, eph.iode, 8);
    assembleFrameBits(out, pos, static_cast<uint32_t>(eph.toe / 300), 11);
    assembleFrameBits(out, pos, getB1CSatType(&eph), 2);
    assembleSignedFrameBits(out, pos, scaleValueToSignedField32(eph.axis - ref_axis, -9), 26);
    assembleSignedFrameBits(out, pos, scaleValueToSignedField32(eph.axis_dot, -21), 25);
    assembleSignedFrameBits(out, pos, scaleValueToSignedField32(eph.delta_n / PI, -44), 17);
    assembleSignedFrameBits(out, pos, scaleValueToSignedField32(eph.delta_n_dot, -57), 23);
    assembleSignedFrameBits(out, pos, scaleValueToSignedField64(eph.m0 / PI, -32), 33);
    assembleFrameBits(out, pos, scaleValueToUnsignedField64(eph.ecc, -34), 33);
    assembleSignedFrameBits(out, pos, scaleValueToSignedField64(eph.w / PI, -32), 33);
}

static void generateEphemeris2Bits(const ephem_t &eph, short out[222])
{
    int pos = 0;

    assembleSignedFrameBits(out, pos, scaleValueToSignedField64(eph.omega0 / PI, -32), 33);
    assembleSignedFrameBits(out, pos, scaleValueToSignedField64(eph.i0 / PI, -32), 33);
    assembleSignedFrameBits(out, pos, scaleValueToSignedField32(eph.omega_dot / PI, -44), 19);
    assembleSignedFrameBits(out, pos, scaleValueToSignedField32(eph.idot / PI, -44), 15);
    assembleSignedFrameBits(out, pos, scaleValueToSignedField32(eph.cis, -30), 16);
    assembleSignedFrameBits(out, pos, scaleValueToSignedField32(eph.cic, -30), 16);
    assembleSignedFrameBits(out, pos, scaleValueToSignedField32(eph.crs, -8), 24);
    assembleSignedFrameBits(out, pos, scaleValueToSignedField32(eph.crc, -8), 24);
    assembleSignedFrameBits(out, pos, scaleValueToSignedField32(eph.cus, -30), 21);
    assembleSignedFrameBits(out, pos, scaleValueToSignedField32(eph.cuc, -30), 21);
}

static void generateClockBits(const ephem_t &eph, short out[69])
{
    int pos = 0;

    assembleFrameBits(out, pos, static_cast<uint32_t>(eph.toc / 300), 11);
    assembleSignedFrameBits(out, pos, scaleValueToSignedField32(eph.af0, -34), 25);
    assembleSignedFrameBits(out, pos, scaleValueToSignedField32(eph.af1, -50), 22);
    assembleSignedFrameBits(out, pos, scaleValueToSignedField32(eph.af2, -66), 11);
}

static int generateSubframe2DataBits(const ephem_t &eph, int bds_week, int how, short data_bits[576])
{
    short eph1[211] = {};
    short eph2[222] = {};
    short clock[69] = {};

    generateEphemeris1Bits(eph, eph1);
    generateEphemeris2Bits(eph, eph2);
    generateClockBits(eph, clock);

    memset(data_bits, 0, B1C_SUBFRAME2_DATA_BITS * sizeof(short));

    int pos = 0;
    assembleFrameBits(data_bits, pos, static_cast<uint32_t>(bds_week), 13);
    assembleFrameBits(data_bits, pos, static_cast<uint32_t>(how), 8);
    assembleFrameBits(data_bits, pos, eph.iodc, 10);
    assembleFrameBitBlock(data_bits, pos, eph1, 211);
    assembleFrameBitBlock(data_bits, pos, eph2, 222);
    assembleFrameBitBlock(data_bits, pos, clock, 69);
    assembleSignedFrameBits(data_bits, pos, scaleValueToSignedField32(eph.tgd_ext[1], -34), 12);
    assembleSignedFrameBits(data_bits, pos, scaleValueToSignedField32(eph.tgd_ext[0] - eph.tgd_ext[1], -34), 12);
    assembleSignedFrameBits(data_bits, pos, scaleValueToSignedField32(eph.tgd_ext[1], -34), 12);
    assembleFrameBits(data_bits, pos, 0, 7);

    return (pos == B1C_SUBFRAME2_DATA_BITS) ? 1 : -1;
}

void b1cSubframe2LdpcEncode600(const short info[600], short coded[1200])
{
    int symbols[B1C_SF2_CODE_SYMBOLS] = {};

    bitsToSymbols6(info, B1C_SF2_INFO_SYMBOLS, symbols);
    ldpcEncodeSymbols(symbols, B1C_SF2_INFO_SYMBOLS, b1c_matrix_gen2);
    symbols6ToBits(symbols, B1C_SF2_CODE_SYMBOLS, coded);
}

static int encodeB1CSubframe2Prepared(bdstime_t g, const ephem_t &prepared_eph, short *subframe)
{
    if (subframe == nullptr || (prepared_eph.valid & 1) == 0)
        return -1;

    memset(subframe, 0, B1C_SUBFRAME2_CODE_BITS_LOCAL * sizeof(short));

    short info[B1C_SUBFRAME2_INFO_BITS_LOCAL] = {};
    const b1c_nav_time_fields_t time_fields = computeB1CNavTimeFields(g);

    if (generateSubframe2DataBits(prepared_eph, g.week, time_fields.how, info) < 0)
        return -2;

    int pos = B1C_SUBFRAME2_DATA_BITS;
    assembleFrameBits(info, pos, crc24qBits(info, B1C_SUBFRAME2_DATA_BITS), 24);
    if (pos != B1C_SUBFRAME2_INFO_BITS_LOCAL)
        return -3;

    b1cSubframe2LdpcEncode600(info, subframe);
    return 1;
}

int generateB1CSubframe2(bdstime_t g, channel_t *chan, ephem_t *eph, short *subframe)
{
    (void)chan;

    if (eph == nullptr)
        return -1;

    ephem_t prepared_eph = *eph;
    prepareB1CEphemeris(prepared_eph);
    return encodeB1CSubframe2Prepared(g, prepared_eph, subframe);
}

static void generateMidiAlmanacBits(const ephem_t &alm, short out[156])
{
    memset(out, 0, 156 * sizeof(short));
    if ((alm.valid & 1) == 0)
        return;

    int pos = 0;
    const double offset = (alm.flag == 1) ? 0.0 : 0.3;

    assembleFrameBits(out, pos, alm.svid, 6);
    assembleFrameBits(out, pos, alm.flag, 2);
    assembleFrameBits(out, pos, static_cast<uint32_t>(alm.week), 13);
    assembleFrameBits(out, pos, static_cast<uint32_t>(alm.toa >> 12), 8);
    assembleFrameBits(out, pos, scaleValueToUnsignedField32(alm.ecc, -16), 11);
    assembleSignedFrameBits(out, pos, scaleValueToSignedField32(alm.i0 / PI - offset, -14), 11);
    assembleFrameBits(out, pos, scaleValueToUnsignedField32(alm.sqrt_a, -4), 17);
    assembleSignedFrameBits(out, pos, scaleValueToSignedField32(alm.omega0 / PI, -15), 16);
    assembleSignedFrameBits(out, pos, scaleValueToSignedField32(alm.omega_dot / PI, -33), 11);
    assembleSignedFrameBits(out, pos, scaleValueToSignedField32(alm.w / PI, -15), 16);
    assembleSignedFrameBits(out, pos, scaleValueToSignedField32(alm.m0 / PI, -15), 16);
    assembleSignedFrameBits(out, pos, scaleValueToSignedField32(alm.af0, -20), 11);
    assembleSignedFrameBits(out, pos, scaleValueToSignedField32(alm.af1, -37), 10);
    assembleFrameBits(out, pos, alm.health, 8);
}

static void generateReducedAlmanacBits(const ephem_t &alm, short out[38])
{
    memset(out, 0, 38 * sizeof(short));
    if ((alm.valid & 1) == 0)
        return;

    int pos = 0;
    const double ref_axis = (alm.flag == 3) ? B1C_MEO_REF_AXIS : B1C_GEO_IGSO_REF_AXIS;

    assembleFrameBits(out, pos, alm.svid, 6);
    assembleFrameBits(out, pos, alm.flag, 2);
    assembleSignedFrameBits(out, pos, scaleValueToSignedField32(alm.sqrt_a * alm.sqrt_a - ref_axis, 9), 8);
    assembleSignedFrameBits(out, pos, scaleValueToSignedField32(alm.omega0 / PI, -6), 7);
    assembleSignedFrameBits(out, pos, scaleValueToSignedField32((alm.m0 + alm.w) / PI, -6), 7);
    assembleFrameBits(out, pos, alm.health, 8);
}

static void generateFallbackAlmanac(const ephem_t *eph, bdstime_t g, ephem_t alm[63])
{
    memset(alm, 0, B1C_ALMANAC_COUNT * sizeof(ephem_t));
    if (eph == nullptr || (eph->valid & 1) == 0 || eph->svid == 0 || eph->svid > 63)
        return;

    ephem_t eph_copy = *eph;
    prepareB1CEphemeris(eph_copy);

    int alm_week = g.week;
    const int tow_seconds = static_cast<int>((static_cast<double>(g.milliseconds) + g.sub_milliseconds) /
                                            BDS_MILLISECONDS_IN_SECOND);
    const int alm_toa = alignBdsAlmanacToa4096(tow_seconds, &alm_week);
    alm[eph_copy.svid - 1] = deriveBdsAlmanacFromEphem(&eph_copy, alm_week, alm_toa);
}

static void generateAlmanacBitTables(
    const ephem_t alm[63],
    int fallback_week,
    int fallback_tow_seconds,
    short midi_bits[63][156],
    short reduced_bits[63][38],
    int &almanac_week,
    int &almanac_toa4096)
{
    almanac_week = fallback_week;
    int aligned_week = fallback_week;
    const int aligned_toa = alignBdsAlmanacToa4096(fallback_tow_seconds, &aligned_week);
    almanac_toa4096 = aligned_toa >> 12;

    for (int i = 0; i < B1C_ALMANAC_COUNT; ++i)
    {
        const ephem_t empty{};
        const ephem_t &item = (alm != nullptr) ? alm[i] : empty;
        generateMidiAlmanacBits(item, midi_bits[i]);
        generateReducedAlmanacBits(item, reduced_bits[i]);

        if ((item.valid & 1) != 0)
        {
            almanac_week = item.week;
            almanac_toa4096 = item.toa >> 12;
        }
    }
}

static int generateSubframe3DataBits(
    int soh,
    uint32_t flags,
    int almanac_week,
    int almanac_toa4096,
    const short midi_bits[63][156],
    const short reduced_bits[63][38],
    short data_bits[240])
{
    memset(data_bits, 0, B1C_SUBFRAME3_DATA_BITS * sizeof(short));

    int pos = 0;
    const int page_index = (soh >= 100) ? (soh - 100) : soh;

    if ((page_index % 10) == 0)
    {
        assembleFrameBits(data_bits, pos, 0, 3);
        assembleFrameBits(data_bits, pos, 1, 3);
        assembleFrameBits(data_bits, pos, 0, 1);
        assembleFrameBits(data_bits, pos, flags, 8);
        assembleFrameBits(data_bits, pos, 0, B1C_SUBFRAME3_DATA_BITS - pos);
    }
    else if ((page_index % 10) == 5)
    {
        assembleFrameBits(data_bits, pos, 0, 3);
        assembleFrameBits(data_bits, pos, 3, 3);
        assembleFrameBits(data_bits, pos, 0, 1);
        assembleFrameBits(data_bits, pos, flags, 8);
        assembleFrameBits(data_bits, pos, 0, B1C_SUBFRAME3_DATA_BITS - pos);
    }
    else if ((page_index % 5) == 1 && page_index < 80)
    {
        const int sv_index = (page_index / 5) * 4;

        assembleFrameBits(data_bits, pos, 0, 3);
        assembleFrameBits(data_bits, pos, 2, 3);
        assembleFrameBits(data_bits, pos, 0, 1);
        assembleFrameBits(data_bits, pos, flags, 8);
        assembleFrameBits(data_bits, pos, 0, 9);
        assembleFrameBits(data_bits, pos, 0, 13);
        assembleFrameBits(data_bits, pos, static_cast<uint32_t>(almanac_week), 13);
        assembleFrameBits(data_bits, pos, static_cast<uint32_t>(almanac_toa4096), 8);

        for (int k = 0; k < 4; ++k)
        {
            if (sv_index + k < B1C_ALMANAC_COUNT)
                assembleFrameBitBlock(data_bits, pos, reduced_bits[sv_index + k], 38);
            else
                assembleFrameBits(data_bits, pos, 0, 38);
        }
    }
    else if (page_index == 99)
    {
        assembleFrameBits(data_bits, pos, 0, 3);
        assembleFrameBits(data_bits, pos, 0, 3);
        assembleFrameBits(data_bits, pos, 0, 1);
        assembleFrameBits(data_bits, pos, flags, 8);
        assembleFrameBits(data_bits, pos, 0, B1C_SUBFRAME3_DATA_BITS - pos);
    }
    else
    {
        const int sv_index = (page_index < 80)
                                 ? (page_index / 5 * 3) + (page_index % 5) - 2
                                 : (page_index / 5 * 4) + (page_index % 5) - 17;

        assembleFrameBits(data_bits, pos, 0, 3);
        assembleFrameBits(data_bits, pos, 4, 3);
        assembleFrameBits(data_bits, pos, 0, 1);
        assembleFrameBits(data_bits, pos, flags, 8);
        assembleFrameBits(data_bits, pos, 0, 9);
        assembleFrameBits(data_bits, pos, 0, 13);

        if (sv_index >= 0 && sv_index < B1C_ALMANAC_COUNT)
            assembleFrameBitBlock(data_bits, pos, midi_bits[sv_index], 156);
        else
            assembleFrameBits(data_bits, pos, 0, 156);
    }

    if (pos < B1C_SUBFRAME3_DATA_BITS)
        assembleFrameBits(data_bits, pos, 0, B1C_SUBFRAME3_DATA_BITS - pos);

    return (pos == B1C_SUBFRAME3_DATA_BITS) ? 1 : -1;
}

static void b1cSubframe3LdpcEncode264(const short info[264], short coded[528])
{
    int symbols[88] = {};

    bitsToSymbols6(info, 44, symbols);
    ldpcEncodeSymbols(symbols, 44, b1c_matrix_gen2);
    symbols6ToBits(symbols, 88, coded);
}

static int encodeB1CSubframe3Prepared(
    bdstime_t g,
    const ephem_t *prepared_eph,
    const ephem_t alm[63],
    short *subframe)
{
    if (subframe == nullptr)
        return -1;

    memset(subframe, 0, B1C_SUBFRAME3_CODE_BITS * sizeof(short));

    short midi_bits[63][156] = {};
    short reduced_bits[63][38] = {};
    short info[B1C_SUBFRAME3_INFO_BITS] = {};

    const int fallback_tow_seconds = static_cast<int>((static_cast<double>(g.milliseconds) + g.sub_milliseconds) /
                                                     BDS_MILLISECONDS_IN_SECOND);
    int almanac_week = 0;
    int almanac_toa4096 = 0;
    generateAlmanacBitTables(alm, g.week, fallback_tow_seconds, midi_bits, reduced_bits,
                          almanac_week, almanac_toa4096);

    const b1c_nav_time_fields_t time_fields = computeB1CNavTimeFields(g);
    if (generateSubframe3DataBits(time_fields.soh,
                               b1cFlagsFromEph(prepared_eph),
                               almanac_week,
                               almanac_toa4096,
                               midi_bits,
                               reduced_bits,
                               info) < 0)
    {
        return -2;
    }

    int pos = B1C_SUBFRAME3_DATA_BITS;
    assembleFrameBits(info, pos, crc24qBits(info, B1C_SUBFRAME3_DATA_BITS), 24);
    if (pos != B1C_SUBFRAME3_INFO_BITS)
        return -3;

    b1cSubframe3LdpcEncode264(info, subframe);
    return 1;
}


int generateB1CSubframe3(bdstime_t g, channel_t *chan, ephem_t *eph, const ephem_t alm[63], short *subframe)
{
    (void)chan;

    ephem_t prepared_eph{};
    ephem_t *prepared_ptr = nullptr;
    if (eph != nullptr)
    {
        prepared_eph = *eph;
        prepareB1CEphemeris(prepared_eph);
        prepared_ptr = &prepared_eph;
    }

    return encodeB1CSubframe3Prepared(g, prepared_ptr, alm, subframe);
}

static void interleaveB1CNavMessage(
    const short subframe1[72],
    const short subframe2[1200],
    const short subframe3[528],
    short nav_msg[1800])
{
    memset(nav_msg, 0, B1C_NAV_BITS * sizeof(short));

    for (int i = 0; i < B1C_SUBFRAME1_BITS; ++i)
        nav_msg[i] = subframe1[i];

    for (int i = 0; i < 11; ++i)
    {
        int p1 = 72 + i * 3;
        int p2 = 73 + i * 3;
        int p3 = 74 + i * 3;

        for (int j = 0; j < 48; ++j)
        {
            nav_msg[p1] = subframe2[i * 96 + j];
            nav_msg[p2] = subframe2[i * 96 + 48 + j];
            nav_msg[p3] = subframe3[i * 48 + j];
            p1 += 36;
            p2 += 36;
            p3 += 36;
        }
    }

    int p1 = 105;
    int p2 = 106;
    int p3 = 107;
    for (int j = 0; j < 48; ++j)
    {
        nav_msg[p1] = subframe2[22 * 48 + j];
        nav_msg[p2] = subframe2[23 * 48 + j];
        nav_msg[p3] = subframe2[24 * 48 + j];
        p1 += 36;
        p2 += 36;
        p3 += 36;
    }
}

int generateBdsB1CMessage(
    bdstime_t g,
    int svid,
    const ephem_t *eph,
    const ephem_t alm[63],
    short nav_msg[1800])
{
    if (eph == nullptr || nav_msg == nullptr || (eph->valid & 1) == 0)
        return -1;

    ephem_t prepared_eph = *eph;
    prepareB1CEphemeris(prepared_eph);

    if (svid <= 0)
        svid = prepared_eph.svid;

    ephem_t fallback[63] = {};
    const ephem_t *almanac = alm;
    if (almanac == nullptr)
    {
        generateFallbackAlmanac(&prepared_eph, g, fallback);
        almanac = fallback;
    }

    short subframe1[B1C_SUBFRAME1_BITS] = {};
    short subframe2[B1C_SUBFRAME2_CODE_BITS_LOCAL] = {};
    short subframe3[B1C_SUBFRAME3_CODE_BITS] = {};

    const b1c_nav_time_fields_t time_fields = computeB1CNavTimeFields(g);
    if (generateSubframe1Bits(svid, time_fields.soh, subframe1) < 0 ||
        encodeB1CSubframe2Prepared(g, prepared_eph, subframe2) < 0 ||
        encodeB1CSubframe3Prepared(g, &prepared_eph, almanac, subframe3) < 0)
    {
        return -1;
    }

    interleaveB1CNavMessage(subframe1, subframe2, subframe3, nav_msg);
    return 1;
}

int generateB1CNavMessage(bdstime_t g, channel_t *chan, ephem_t *eph, const ephem_t alm[63])
{
    if (chan == nullptr || eph == nullptr || chan->nav_bit == nullptr)
        return -1;

    short nav_msg[B1C_NAV_BITS] = {};
    if (generateBdsB1CMessage(g, eph->svid, eph, alm, nav_msg) < 0)
        return -1;

    for (int i = 0; i < B1C_NAV_BITS; ++i)
        chan->nav_bit[i] = nav_msg[i] ? static_cast<short>(-1) : static_cast<short>(1);

    return 1;
}

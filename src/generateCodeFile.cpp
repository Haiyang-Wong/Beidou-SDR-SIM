/**
 * \file generateCodeFile.cpp
 * \brief B1C 主码生成相关函数实现
 * \author LackWood Du
 * \date 2025-12-19
 */

#include "bds_sdr.h"
#pragma message("Compiling generateCodeFile.cpp")

/**
 * \brief 计算模幂
 * \note Euler 判别，判断base是否为mod下的二次剩余
 */
long long modPow(long long base, long long exp, long long mod)
{
    long long result = 1;
    base %= mod;

    while (exp > 0)
    {
        if (exp & 1)
            result = (result * base) % mod;
        base = (base * base) % mod;
        exp >>= 1;
    }
    return result;
}

/**
 * \brief 计算勒让德符号
 * \note 计算 a 关于素数 p 的勒让德符号
 */
int legendreSymbol(int a, int p)
{
    if (a == 0)
        return 0;

    long long t = modPow(a, (p - 1) / 2, p);

    if (t == 1)
        return 1;
    else if (t == p - 1)
        return 0;
    // -1;
    else
        throw std::runtime_error("Unexpected value in Legendre symbol");
}

/**
 * \brief 判断整数是否为素数
 */
bool isPrimeInt(int n)
{
    if (n <= 1)
        return false;
    if (n == 2)
        return true;
    if (n % 2 == 0)
        return false;
    for (int d = 3; 1LL * d * d <= n; d += 2)
    {
        if (n % d == 0)
            return false;
    }
    return true;
}

/**
 * \brief 生成勒让德序列
 * \param[in] p 奇素数， 默认值为10243
 * \note 生成长度为p的勒让德序列，p为奇素数
 * \return 生成的勒让德序列
 */
std::vector<int> generateLegendreSequence(int p = 10243)
{
    if ((p % 2 == 0) || !isPrimeInt(p))
    {
        throw std::invalid_argument("p must be an odd prime");
    }

    std::vector<int> seq(p);
    for (int k = 0; k < p; ++k)
    {
        seq[k] = legendreSymbol(k, p);
    }
    return seq;
}

/**
 * \brief 生成Weil序列
 * \param[in] L 勒让德序列
 * \param[in] w 移位值
 * \note 根据给定的勒让德序列L和移位值w，生成Weil序列W
 * \return 生成的Weil序列
 */
std::vector<int> generateWeilSequence(
    const std::vector<int> &L,
    int w)
{
    int N = L.size();
    std::vector<int> W(N);

    for (int k = 0; k < N; ++k)
    {
        int kp = (k + w) % N;
        W[k] = L[k] ^ L[kp]; // XOR
    }
    return W;
}

/**
 * \brief 按 ICD 中的定义截取主码
 * \param[in] W Weil序列，不同的卫星对应不同的Weil序列
 * \param[in] p 截取点，取值范围是 1 - N（1-based）
 * \param[in] codeLength 码片长度，默认值为10230
 * \return 生成的主码序列
 */
std::vector<int> generateB1CPrimaryCode(
    const std::vector<int> &W,
    int p, // 1-based
    int codeLength)
{
    int N = W.size();
    std::vector<int> code(codeLength);

    for (int n = 0; n < codeLength; ++n)
    {
        int index = (n + p - 1) % N; // 注意 p-1
        code[n] = W[index];
    }
    return code;
}

/**
 * \brief 将码片序列转换为字符串形式
 * \param[in] code 码片序列
 * \return 转换后的字符串
 */
std::string codeToString(const std::vector<int> &code)
{
    std::string s;
    s.reserve(code.size());

    for (int c : code)
    {
        s.push_back(c ? '1' : '0');
    }
    return s;
}

/**
 * \brief 将生成的码片序列写入CSV文件
 * \param[in] file 输出文件流
 * \param[in] prn 卫星编号
 * \param[in] code 码片序列
 */
void writeOneCodeRow(
    std::ofstream &file,
    int prn,
    const std::vector<int> &code)
{
    file << prn << "," << codeToString(code) << "\n";
}

/**
 * \brief 将24个码片转换为8进制字符串
 * \param[in] chips 长度为24的码片数组
 * \return 转换后的8进制字符串
 */
std::string chips24ToOctal(const std::vector<int> &chips)
{
    if (chips.size() != 24)
    {
        throw std::runtime_error("chips24ToOctal: input size must be 24");
    }

    std::ostringstream oss;

    for (int i = 0; i < 24; i += 3)
    {
        int value = (chips[i] << 2) | (chips[i + 1] << 1) | (chips[i + 2]);
        oss << value;
    }

    return oss.str();
}

/**
 * \brief 验证生成的主码是否符合ICD标准
 * \param[in] prn 卫星编号
 * \param[in] code 生成的主码序列
 * \param[in] refHead 参考头部8进制字符串
 * \param[in] refTail 参考尾部8进制字符串
 */
void verifyAgainstICD(
    int prn,
    const std::vector<int> &code,
    const std::string &refHead,
    const std::string &refTail)
{
    const int L = code.size();

    std::vector<int> head(code.begin(), code.begin() + 24);
    std::vector<int> tail(code.end() - 24, code.end());

    std::string head_oct = chips24ToOctal(head);
    std::string tail_oct = chips24ToOctal(tail);

    if (head_oct != refHead || tail_oct != refTail)
    {
        std::ostringstream err;
        err << "B1C primary code verification FAILED for PRN "
            << prn << "\n"
            << "Expected head: " << refHead << ", got: " << head_oct << "\n"
            << "Expected tail: " << refTail << ", got: " << tail_oct;

        throw std::runtime_error(err.str());
    }
}

/**
 * \brief 主生成函数
 * \param[in] params 包含卫星参数的结构体数组
 * \param[in] legendre 勒让德序列
 * \param[in] outputPrefix 输出文件名前缀
 * \param[in] codeLength 码片长度，默认值为10230
 */
void generateB1CComponentCodes(
    const std::vector<B1CParams> &params,
    const std::vector<int> &legendre,
    const std::string &outputCsvPath,
    const int codeLength // weil序列长度
)
{
    std::ofstream file(outputCsvPath);
    if (!file)
    {
        throw std::runtime_error("Cannot open output CSV: " + outputCsvPath);
    }

    // CSV 表头
    file << "PRN,Code\n";

    for (const auto &item : params)
    {

        // 1. Weil 码
        auto W = generateWeilSequence(legendre, item.w);

        // 2. 截断生成主码，指定主码 / 主码的长度
        auto code = generateB1CPrimaryCode(W, item.p, codeLength);

        // 3. ICD 校验
        verifyAgainstICD(
            item.prn,
            code,
            item.head24_oct,
            item.tail24_oct);

        // 4. 写入一行
        writeOneCodeRow(file, item.prn, code);
    }
}

/**
 * \brief 从CSV文件加载B1C参数
 * \param[in] csvPath CSV文件路径
 * \return 加载的B1C参数数组
 */
std::vector<B1CParams> loadB1CParamsFromCSV(const std::string &csvPath)
{
    std::ifstream file(csvPath);
    if (!file.is_open())
    {
        throw std::runtime_error("Cannot open CSV file: " + csvPath);
    }

    std::vector<B1CParams> params;
    std::string line;

    // 1. 读取并丢弃表头
    if (!std::getline(file, line))
    {
        throw std::runtime_error("CSV file is empty: " + csvPath);
    }

    // 2. 逐行解析数据
    int lineNumber = 1;
    while (std::getline(file, line))
    {
        ++lineNumber;

        if (line.empty())
        {
            continue; // 跳过空行
        }

        std::stringstream ss(line);
        std::string field;
        std::vector<std::string> fields;

        while (std::getline(ss, field, ','))
        {
            fields.push_back(field);
        }

        if (fields.size() != 5)
        {
            std::ostringstream err;
            err << "Invalid CSV format at line " << lineNumber
                << " (expected 5 fields, got " << fields.size() << ")";
            throw std::runtime_error(err.str());
        }

        B1CParams p;
        try
        {
            p.prn = std::stoi(fields[0]);
            p.w = std::stoi(fields[1]);
            p.p = std::stoi(fields[2]);
            p.head24_oct = fields[3];
            p.tail24_oct = fields[4];
        }
        catch (const std::exception &)
        {
            std::ostringstream err;
            err << "Failed to parse numeric fields at line " << lineNumber;
            throw std::runtime_error(err.str());
        }

        params.push_back(p);
    }

    if (params.empty())
    {
        throw std::runtime_error("No valid B1C parameters found in CSV: " + csvPath);
    }

    return params;
}

/**
 * \brief 修正CSV文件中第4、5列的8进制字段宽度，不足8位时前面补零
 * \param[in] csvPath CSV文件路径
 */
void fixOctalFieldWidthInCSV(const std::string &csvPath)
{
    std::ifstream inFile(csvPath);
    if (!inFile.is_open())
    {
        throw std::runtime_error("Cannot open CSV file for reading: " + csvPath);
    }

    std::vector<std::string> lines;
    std::string line;

    // 1. 读入所有行
    while (std::getline(inFile, line))
    {
        lines.push_back(line);
    }
    inFile.close();

    if (lines.empty())
    {
        return; // 空文件直接返回
    }

    // 2. 从第 2 行开始处理（跳过表头）
    for (size_t i = 1; i < lines.size(); ++i)
    {
        if (lines[i].empty())
        {
            continue;
        }

        std::stringstream ss(lines[i]);
        std::vector<std::string> fields;
        std::string field;

        while (std::getline(ss, field, ','))
        {
            fields.push_back(field);
        }

        // 必须至少有 5 列
        if (fields.size() < 5)
        {
            continue; // 或选择抛异常
        }

        // 第 4、5 列补零（索引 3,4）
        for (int idx : {3, 4})
        {
            if (fields[idx].length() < 8)
            {
                fields[idx].insert(
                    fields[idx].begin(),
                    8 - fields[idx].length(),
                    '0');
            }
        }

        // 重新拼接该行
        std::ostringstream rebuilt;
        for (size_t j = 0; j < fields.size(); ++j)
        {
            if (j > 0)
                rebuilt << ",";
            rebuilt << fields[j];
        }

        lines[i] = rebuilt.str();
    }

    // 3. 写回原文件
    std::ofstream outFile(csvPath, std::ios::trunc);
    if (!outFile.is_open())
    {
        throw std::runtime_error("Cannot open CSV file for writing: " + csvPath);
    }

    for (const auto &l : lines)
    {
        outFile << l << "\n";
    }
}
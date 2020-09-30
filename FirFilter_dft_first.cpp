#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <charconv>
#include <complex>

constexpr int TRANSFORM_POINTS = 10000;

namespace user_literal
{
    constexpr double pi = 3.14159265358979323846264338327950288;
}

constexpr std::complex<double> operator"" _i_comp(const long double val)
{
    std::complex<double> complex_val(0.0, static_cast<double>(val));
    return complex_val;
}

constexpr std::complex<double> operator"" _r_comp(const long double val)
{
    std::complex<double> complex_val(static_cast<double>(val), 0.0);
    return complex_val;
}

const std::string output_file_name = "out.txt";

std::vector<std::complex<double>> convolve_freq(const std::vector<std::complex<double>>& vec1, const std::vector<std::complex<double>>& vec2);
std::vector<double> convolve(const std::vector<double>& vec1, const std::vector<double>& vec2);
std::vector<std::complex<double>> dft(const std::vector<double>& data_vec, const int transform_points);

/*
 * 0 : name of file
 * 1 : signal file
 * 2 : filter file 
 * 3 : output file name
 */
int main(int argv, char** argc)
{
    if(argv != 4 && argv != 2)
    {
        std::cout << "File name arguments are invalid.Please run this program with --help command to confirm arguments." << std::endl;
        return -1;
    }
    const std::string arg1_string = *(argc+1);
    if(argv == 2 && arg1_string.find("--help") != std::string::npos)
    {
        std::cout << "Usage: .\\a.exe [signal file path] [filter file path] [output file name]"  << std::endl;
        return 0;
    }

    const char* const sig_file_name = *(argc + 1);
    const char* const filter_file_name = *(argc + 2);
    const char* const output_file_name = *(argc + 3);
    std::ifstream sig_file, filter_file;

    /*open files*/
    sig_file.open(sig_file_name);
    filter_file.open(filter_file_name);

    if(!sig_file)
    {
        std::cout << "\"" << sig_file_name << "\"" << " is can not be opened." << std::endl;
    }
    if(!filter_file)
    {
        std::cout << "\"" << filter_file_name << "\"" << " is can not be opened." << std::endl;
    }

    std::vector<double> sig_data_vec = {};
    std::vector<double> filter_data_vec = {};

    std::string str_buf;
    double val = 0.0;
    while(std::getline(sig_file, str_buf))
    {
        sig_data_vec.push_back(std::stod(str_buf));
    }
    while(std::getline(filter_file, str_buf))
    {
        filter_data_vec.push_back(std::stod(str_buf));
    }

    auto&& transformed_filter = dft(filter_data_vec, TRANSFORM_POINTS);
    auto&& transformed_sig = dft(sig_data_vec, TRANSFORM_POINTS);

    auto&& convolved_data = convolve_freq(transformed_sig, transformed_filter);

    std::ofstream output_file(output_file_name, std::ios::trunc);

    for(auto& i : convolved_data)
    {
        output_file << std::to_string(20.0 * std::log10(std::abs(i))) << std::endl;
    }

    return 0;
}

std::vector<double> convolve(const std::vector<double>& vec1, const std::vector<double>& vec2)
{
    std::vector<double> convolved_vec(vec1.size() + vec2.size() - 1);
    int i = 0;
    for(auto& convolved_vec_element : convolved_vec)
    {
        convolved_vec_element = 0;
        int initial_p_val = (1+i-vec1.size());
        if(initial_p_val < 0)initial_p_val = 0; 
        for(int p = initial_p_val; p <= i && vec2.size() > p; ++p)
        {
            if((i-p) >= 0)convolved_vec_element += vec1.at(i-p) * vec2.at(p);
        }
        ++i;
    }
    return convolved_vec;
}

std::vector<std::complex<double>> convolve_freq(const std::vector<std::complex<double>>& vec1, const std::vector<std::complex<double>>& vec2)
{
    std::vector<std::complex<double>> convolved_vec(vec1.size() > vec2.size() ? vec1.size() : vec2.size());
    int count = 0;
    for(auto& i : convolved_vec)
    {
        if(count >= vec1.size() || count >= vec2.size())
        {
            i = 0;
        }
        else
        {
            i = vec1.at(count) * vec2.at(count);
        }
        ++count;
    }
    return convolved_vec;
}

std::vector<std::complex<double>> dft(const std::vector<double>& data_vec, const int transform_points)
{
    std::vector<std::complex<double>> transformed_vec(transform_points);
    int count = 0;
    const double inv_N = 1.0 / static_cast<double>(transform_points);
    for(auto& ele : transformed_vec)
    {
        ele = {};
        for(int i = 0; i < transform_points && i < data_vec.size(); ++i)
        {
            ele += data_vec.at(i) * (1.0_r_comp*std::cos(-2.0*user_literal::pi*inv_N*static_cast<double>(count*i)) + 1.0_i_comp*std::sin(-2.0*user_literal::pi*inv_N*static_cast<double>(count*i)) );
        }
        ele *= inv_N;
        ++count;
    }
    return transformed_vec;
}
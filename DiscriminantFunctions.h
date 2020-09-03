#ifndef DISCRIMINANTFUNCTIONS_H
#define DISCRIMINANTFUNCTIONS_H

#include <vector>
#include <unordered_map>
#include <string>
#include <map>

std::vector<long double> datafile_reader(const std::string & filename);
long double discriminant_function(long double a, long double b, long double M2, long double M4, long double size_n);
std::vector<long double> initial_moment_finder(long double input_data [], int size);
std::vector<long double> initial_moment_finder(const std::vector<long double> & input_data, int size);
std::vector<long double> moment_updater(long double old_val, long double new_val, int size, std::vector<long double> initial_moments);
std::vector<std::string> split(const std::string & s, char c);
std::unordered_map<std::string, int> header_parser(const std::string & header, char c);
std::map<std::string, std::vector<long double>> vs_file_reader(const char* filename);
std::vector<long double> discriminant_values_generator(const std::vector<long double>& input_data, int window_size);
std::unordered_map<std::string, std::string> config_reader(const char* filename);
#endif

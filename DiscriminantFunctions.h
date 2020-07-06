#ifndef DISCRIMINANTFUNCTIONS_H
#define DISCRIMINANTFUNCTIONS_H

#include <vector>
#include <unordered_map>
#include <string>
#include <map>

std::vector<double> datafile_reader(const char* filename);
double discriminant_function(double a, double b, double M2, double M4, double size_n);
std::vector<double> initial_moment_finder(double input_data [], int size);
std::vector<double> initial_moment_finder(const std::vector<double> & input_data, int size);
std::vector<double> moment_updater(double old_val, double new_val, int size, std::vector<double> initial_moments);
std::vector<std::string> split(const std::string & s, char c);
std::unordered_map<std::string, int> header_parser(const std::string & header, char c);
std::map<std::string, std::vector<double>> vs_file_reader(const char* filename);
std::vector<double> discriminant_values_generator(const std::vector<double>& input_data, int window_size);

#endif

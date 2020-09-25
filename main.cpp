#include <iostream>
#include <iterator>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <unordered_map>
#include <map>
#include <unordered_set>
#include <typeinfo>
#include <cmath>
#include <chrono>
#include "DiscriminantFunctions.h"
#include "FilterFunctions.h"

using namespace std;
using namespace std::chrono;
void runner(const string & filename, const string &  output_file_name, const string &  win_size_discrim, const string &  win_size_SG, const string &  polynom, const string &  cutoff, const string & tolerance, const string &  secondary_file_path, const string & auto_window){

	int poly_deg = std::stoi(polynom);
	long double mosaic_threshold = std::stold(cutoff);
	long double tol_val = -1.0L * std::stold(tolerance);
	int discrim_win_size;
	int sg_win_size;
	std::vector<long double> input_data = datafile_reader(filename);
	if (auto_window.compare("Y") == 0){
		// auto window is set to true
		int input_data_size = input_data.size();
		long double base_win_size = floor((long double)input_data_size / 100.0L);
		int base_int_size = (int)base_win_size;
		if (base_int_size % 2 == 0){
			discrim_win_size = base_int_size;
			sg_win_size = base_int_size + 1;
		}
		else{
			discrim_win_size = base_int_size - 1;
			sg_win_size = base_int_size;
		}
	}
	else{
		// convert the two window sizes into integers
		discrim_win_size = std::stoi(win_size_discrim);
		sg_win_size = std::stoi(win_size_SG);
	}
	// take the first x number of values and calculate initital moments
	std::cout << "Discriminant window size is " << discrim_win_size << std::endl;
	std::cout << "SG window size is " << sg_win_size << std::endl;

	int slide_iterations = input_data.size() - discrim_win_size;
	if (slide_iterations < 1){
		std::cout << "window size of " << discrim_win_size << " is bigger than the size of data " << input_data.size() << std::endl;
	}
	else{
		std::vector<long double> initial_datapoints(discrim_win_size, 0.0L);
		for (int i = 0; i != discrim_win_size; ++i){
			initial_datapoints.at(i) = input_data.at(i);
		}
		std::vector<long double> moments = initial_moment_finder(initial_datapoints, discrim_win_size);

		// calculate the first discriminant value and store it into a vector of discriminant values
		std::vector<long double> discriminants(slide_iterations + 1, 0.0L);
		discriminants.at(0) = discriminant_function(200.0L, 0.05L, moments.at(1), moments.at(3), discrim_win_size);

		// calculate the rest of the discriminants using a sliding window
		int old_value_index = 0;
		int new_value_index = discrim_win_size;
		for (int i = 0; i != slide_iterations; ++i){
			// first calculate new moments
			moments = moment_updater(input_data[old_value_index], input_data[new_value_index], discrim_win_size, moments);
			// calculate the new discrimanant and add to vector
			discriminants.at(i + 1) = discriminant_function(200.0L, 0.05L, moments.at(1), moments.at(3), discrim_win_size);
			++old_value_index;
			++new_value_index;
		}

		// apply SG filter to data
		std::vector<long double> discrim_filtered = sg_smooth(discriminants, sg_win_size, poly_deg);

		// write the output filtered values to file
		std::ofstream result_file(output_file_name);
		std::ofstream secondary_file(secondary_file_path);
		int counter = 1;
		int counter2 = discrim_win_size;
		bool isMosaic = true;
		long double diff = 0.0L;
		bool mosaicismPresent = false;
		if (result_file.is_open() && secondary_file.is_open()){
			// write out columns headers
			result_file << "Start_Pos" << "\t" << "End_Pos"<< "\t" << "Discriminant_Value" << "\t" << "Is_Mosaic" << "\n";
			secondary_file << "Start_Pos" << "\t" << "End_Pos"<< "\t" << "Discriminant_Value" << "\t" << "Distance from cutoff" << "\n";
			for (auto & elem: discrim_filtered){
				isMosaic = (elem > mosaic_threshold) ? true : false;
				result_file << counter << "\t" <<  counter2 << "\t" << elem << "\t" << isMosaic << "\n";
				diff = elem - mosaic_threshold;
				if (isMosaic){
					secondary_file << counter << "\t" <<  counter2 << "\t" << elem << "\t" << diff << "\n";
					if (!mosaicismPresent){
						mosaicismPresent = true;
					}
				}
				else{
					if (diff > tol_val){
						secondary_file << counter << "\t" <<  counter2 << "\t" << elem << "\t" << diff << "\n";
					}
				}
				++counter;
				++counter2;
			}
			result_file.close();
			secondary_file.close();
		}
		else {
			std::cout << "Unable to write to file" << std::endl;
		}
		if (mosaicismPresent){
			cout << "Mosaicism was detected" << endl;
		}
		else {
			cout << "Mosaicism was not detected" << endl;
		}
	}

}

int main(int argc, char* argv[]){
	// arg 1: input_file path, arg 2: output_file path, arg 3: discriminant window size, arg 4: width of SG filter, arg 5: polynomial degree of SG filter,
	// arg 6: threshold for demeaning an area mosaic, arg 7: tolerance value, arg 8: secondary file name

	std::unordered_map<std::string, std::string> config_vals = config_reader(argv[1]);
	std::string arg1 = config_vals.at("INPUT_FILE_PATH");
	std::string arg2 = config_vals.at("OUTPUT_FILE_PATH");
	std::string arg3 = config_vals.at("DISCRIMINANT_WINDOW_SIZE");
	std::string arg4 = config_vals.at("SG_FILTER_WIDTH");
	std::string arg5 = config_vals.at("SG_FILTER_POLYNOMIAL_DEGREE");
	std::string arg6 = config_vals.at("MOSAIC_THRESHOLD");
	std::string arg7 = config_vals.at("TOLERANCE_VALUE");
	std::string arg8 = config_vals.at("SECONDARY_FILE_PATH");
	std::string arg9 = config_vals.at("AUTO_WINDOW_SIZE");
	auto start = high_resolution_clock::now();
	runner(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9);
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(stop - start);
	cout << "Time taken by function: " << duration.count() << " milliseconds" << endl;
	return 0;

}

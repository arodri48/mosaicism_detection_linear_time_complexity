#include <iostream>
#include <iterator>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <unordered_map>
#include <map>
#include <unordered_set>

#include "DiscriminantFunctions.h"
#include "FilterFunctions.h"

void runner(const char* filename, const char* output_file_name, const char* win_size_discrim, const char* win_size_SG, const char* polynom, const char* cutoff){
	// convert the two window sizes into integers
	int discrim_win_size = std::stoi(win_size_discrim);
	int sg_win_size = std::stoi(win_size_SG);
	int poly_deg = std::stoi(polynom);
	double mosaic_threshold = std::stod(cutoff);

	// read in the input data into memory
	std::vector<double> input_data = datafile_reader(filename);

	// take the first x number of values and calculate initital moments
	int slide_iterations = input_data.size() - discrim_win_size;
	std::vector<double> initial_datapoints(discrim_win_size, 0.0);
	for (int i = 0; i != discrim_win_size; ++i){
		initial_datapoints.at(i) = input_data.at(i);
	}
	std::vector<double> moments = initial_moment_finder(initial_datapoints, discrim_win_size);

	// calculate the first discriminant value and store it into a vector of discriminant values
	std::vector<double> discriminants(slide_iterations + 1, 0.0);
	discriminants.at(0) = discriminant_function(200.0, 0.05, moments.at(1), moments.at(3), discrim_win_size);

	// calculate the rest of the discriminants using a sliding window
	int old_value_index = 0;
	int new_value_index = discrim_win_size;
	for (int i = 0; i != slide_iterations; ++i){
		// first calculate new moments
		moments = moment_updater(input_data[old_value_index], input_data[new_value_index], discrim_win_size, moments);
		// calculate the new discrimanant and add to vector
		discriminants.at(i + 1) = discriminant_function(200.0, 0.05, moments.at(1), moments.at(3), discrim_win_size);
		++old_value_index;
		++new_value_index;
	}

	// apply SG filter to data
	std::vector<double> discrim_filtered = sg_smooth(discriminants, sg_win_size, poly_deg);

	// write the output filtered values to file
	std::ofstream result_file(output_file_name);
	int counter = 1;
	int counter2 = discrim_win_size;
	bool isMosaic = true;
	if (result_file.is_open()){
		// write out columns headers
		result_file << "Start_Pos" << "\t" << "End_Pos"<< "\t" << "Discriminant_Value" << "\t" << "Is_Mosaic" << "\n";
		for (auto & elem: discrim_filtered){
			isMosaic = (elem > mosaic_threshold) ? true : false;
			result_file << counter << "\t" <<  counter2 << "\t" << elem << "\t" << isMosaic << "\n";
			++counter;
			++counter2;
		}
		result_file.close();
	}
	else {
		std::cout << "Unable to write to file";
	}
}

int main(int argc, char* argv[]){
	// arg 1: input_file path, arg 2: output_file path, arg 3: discriminant window size, arg 4: width of SG filter, arg 5: polynomial degree of SG filter,
	// arg 6: threshold for demeaning an area mosaic
	runner(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6]);
	return 0;
}

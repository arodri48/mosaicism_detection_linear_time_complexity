#include <iostream>
#include <iterator>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <unordered_map>
#include <map>
#include <unordered_set>

std::vector<long double> datafile_reader(const char* filename){
	// assumes a input file with a single column of data with no headers
	std::ifstream is(filename);
	std::istream_iterator<long double> start(is), end;
	std::vector<long double> numbers(start, end);
	return numbers;
}
long double discriminant_function(long double a, long double b, long double M2, long double M4, long double size_n){
	long double nm1 = size_n - 1.0L;
	long double m2_sqr_n = M2 * M2 * size_n;
	long double discrim_val = ((a* m2_sqr_n * M2) - (b * M4 * nm1 * nm1 * nm1))/(nm1*m2_sqr_n) + 3.0L*b;
	return discrim_val;
}
std::vector<long double> initial_moment_finder(long double input_data [], int size){
	std::vector<long double> moments(4,0.0L);
	int num_elem = size;
	long double n = 0.0L;
	long double inv_n = 0.0L;
	long double val = 0.0L;
	long double delta = 0.0L;
	long double A = 0.0L;
	long double B = 0.0L;
	long double mean = 0.0L;
	long double mom2 = 0.0L;
	long double mom3 = 0.0L;
	long double mom4 = 0.0L;
	for (int i = 0; i != num_elem; ++i){
		n += 1.0L;
		inv_n = 1.0L / n;
		val = input_data[i];
		delta = val - mean;
		A = delta * inv_n;
		mean += A;
		mom4 += A * (A * A * delta * i * (n * (n - 3.0L) + 3.0L) + (6.0L * A * mom2) - (4.0L * mom3));
		B = val - mean;
    mom3 += A * (B * delta * (n - 2.0L) - 3.0L * mom2);
    mom2 += delta * B;
	}
	moments.at(0) = mean;
	moments.at(1) = mom2;
	moments.at(2) = mom3;
	moments.at(3) = mom4;
	return moments;
}
std::vector<long double> initial_moment_finder(const std::vector<long double> & input_data, int size){
	std::vector<long double> moments(4,0.0L);
	int num_elem = size;
	long double n = 0.0L;
	long double inv_n = 0.0L;
	long double val = 0.0L;
	long double delta = 0.0L;
	long double A = 0.0L;
	long double B = 0.0L;
	long double mean = 0.0L;
	long double mom2 = 0.0L;
	long double mom3 = 0.0L;
	long double mom4 = 0.0L;
	for (int i = 0; i != num_elem; ++i){
		n += 1.0L;
		inv_n = 1.0L / n;
		val = input_data.at(i);
		delta = val - mean;
		A = delta * inv_n;
		mean += A;
		mom4 += A * (A * A * delta * i * (n * (n - 3.0L) + 3.0L) + 6.0L * A * mom2 - 4.0L * mom3);
		B = val - mean;
    mom3 += A * (B * delta * (n - 2.0L) - 3.0L * mom2);
    mom2 += delta * B;
	}
	moments.at(0) = mean;
	moments.at(1) = mom2;
	moments.at(2) = mom3;
	moments.at(3) = mom4;
	return moments;
}
std::vector<long double> moment_updater(long double old_val, long double new_val, int size, std::vector<long double> initial_moments){
	long double nm1 = (long double)(size) - 1.0L;
	long double nm2 = (long double)(size) - 2.0L;
	long double old_val_min_mu = old_val - initial_moments[0];
	long double mom2 = initial_moments[1];
	long double mom3 = initial_moments[2];
	long double mom4 = initial_moments[3];
	long double mean = 0.0L;
	long double nm1_sqr = nm1 * nm1;
	long double old_val_min_mu_sqr = old_val_min_mu * old_val_min_mu;
	// calculate 2nd moment with old value removed
	mom2 -= (size/nm1) * old_val_min_mu_sqr ;
	// calculate 3rd moment with old value removed
	mom3 -= (nm2 * size * old_val_min_mu_sqr - 3.0L * mom2 * nm1) * (old_val_min_mu / nm1_sqr);
	// calculate 4th moment with old value removed
	mom4 -= (old_val_min_mu / (nm1_sqr * nm1))* ((6.0L * mom2 * nm1 * old_val_min_mu) - (4.0L * mom3 * nm1_sqr) + (size * (size * size - 3.0L * size + 3.0L)*old_val_min_mu_sqr * old_val_min_mu));
	// calculate the mean with old value removed
	mean = (size * initial_moments[0] - old_val) / nm1;
	// update the moments with the new value
	long double delta = new_val - mean;
	long double inv_n = 1.0L / (long double)(size);
	long double A = delta * inv_n;
	mean += A;
	mom4 += A * (A * A * delta * nm1 * (size * (size - 3.0L) + 3.0L) + 6.0L * A * mom2 - 4.0L * mom3);
	long double B = new_val - mean;
	mom3 += A * (B * delta * nm2 - 3.0L * mom2);
	mom2 += delta * B;
	std::vector<long double> new_moments(4,0.0L);
	new_moments.at(0) = mean;
	new_moments.at(1) = mom2;
	new_moments.at(2) = mom3;
	new_moments.at(3) = mom4;
	return new_moments;
}
std::vector<std::string> split(const std::string & s, char c){
	std::vector<std::string> splitted;
	std::string word;
	for (char ch : s){
		if ((ch == c) && (!word.empty())){
			splitted.push_back(word);
			word.clear();
		}
		else{
			word += ch;
		}
	}
	if(!word.empty()){
		splitted.push_back(word);
	}
	return splitted;
}
std::unordered_map<std::string, int> header_parser(const std::string & header, char c){
	std::vector<std::string> header_split = split(header, c);
	int size = header_split.size();
	std::unordered_map<std::string, int> header_indexes;
	for (int i = 0; i != size; ++i){
		header_indexes.insert(std::pair<std::string, int>(header_split.at(i), i));
	}
	return header_indexes;
}
std::map<std::string, std::vector<long double>> vs_file_reader(const char* filename){
	//open file up for reading
	std::string file_path(filename);
	std::ifstream in_s(file_path);
	std::string line;
	std::vector<std::string> line_split;
	std::map<std::string, std::vector<long double>> chr_data;
	// read the header in and parse the header
	if (in_s){
		getline(in_s, line);
		std::unordered_map<std::string, int> header_indexes = header_parser(line, '\t'); // O(1) access
		// read each line, split it, check if chromosome name already exists in a set. If so, add to current vector. If new chromosome name, push current vector into map, clear the vector, add the new name to the set,
		// and add new value
		int chr_index = header_indexes.at("Chr");
		int data_index = header_indexes.at("Data");
		std::unordered_set<std::string> chr_names;
		std::vector<long double> cur_chr_vals;
		std::string old_chr_name;
		while(getline(in_s, line)){
			line_split = split(line, '\t'); // split line
			// check if chromosome name already exists
			if (chr_names.find(line_split.at(chr_index)) != chr_names.end()){// this is case if chromosome already exists
				cur_chr_vals.push_back(std::stod(line_split.at(data_index)));
			}
			else{ // this is case if chromosome is new name
				// check if set is empty
				if (!chr_names.empty()){// in the case that this is not the first chromosome name
					chr_names.insert(line_split.at(chr_index)); // insert name into chr_names
					chr_data.insert(std::pair<std::string, std::vector<long double>> (old_chr_name, cur_chr_vals));
					cur_chr_vals.clear();
					cur_chr_vals.push_back(std::stod(line_split.at(data_index)));
					old_chr_name = line_split.at(chr_index);
				}
				else{
					old_chr_name = line_split.at(chr_index);
					chr_names.insert(line_split.at(chr_index)); // insert name into chr_names
					cur_chr_vals.push_back(std::stod(line_split.at(data_index)));
				}
			}
		}
		// do one final insertion into chr_data after reading the final line
		chr_data.insert(std::pair<std::string, std::vector<long double>> (old_chr_name, cur_chr_vals));
		in_s.close();
	}
	else {
		std::cout << "Could not open: " << file_path << std::endl;
	}
	return chr_data;
}
std::vector<long double> discriminant_values_generator(const std::vector<long double> & input_data, int window_size){
	int slide_iterations = input_data.size() - window_size;
	std::vector<long double> initial_datapoints(window_size, 0.0L);
	for (int i = 0; i != window_size; ++i){
		initial_datapoints.at(i) = input_data.at(i);
	}
	std::vector<long double> moments = initial_moment_finder(initial_datapoints, window_size);


	// calculate the first discrimanant value and store it into a vector of discrimanant values
	std::vector<long double> discriminants(slide_iterations + 1, 0.0L);
	discriminants.at(0) = discriminant_function(200.0L, 0.05L, moments.at(1), moments.at(3), window_size);


	// implement the sliding window
	int old_value_index = 0;
	int new_value_index = window_size;
	//std::vector<long double> new_moments;
	for (int i = 0; i != slide_iterations; ++i){
		// first calculate new moments
		moments = moment_updater(input_data[old_value_index], input_data[new_value_index], window_size, moments);
		// overwrite the old moments with the new ones
		//moments = new_moments;
		// calculate the new discrimanant and add to vector
		discriminants.at(i + 1) = discriminant_function(200.0L, 0.05L, moments.at(1), moments.at(3), window_size);
		++old_value_index;
		++new_value_index;
	}
	return discriminants;
}

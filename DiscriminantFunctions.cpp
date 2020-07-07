#include <iostream>
#include <iterator>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <unordered_map>
#include <map>
#include <unordered_set>

std::vector<double> datafile_reader(const char* filename){
	// assumes a input file with a single column of data with no headers
	std::ifstream is(filename);
	std::istream_iterator<double> start(is), end;
	std::vector<double> numbers(start, end);
	return numbers;
}
double discriminant_function(double a, double b, double M2, double M4, double size_n){
	double nm1 = size_n - 1.0;
	double m2_sqr_n = M2 * M2 * size_n;
	double discrim_val = ((a* m2_sqr_n * M2) - (b * M4 * nm1 * nm1 * nm1))/(nm1*m2_sqr_n) + 3.0*b;
	return discrim_val;
}
std::vector<double> initial_moment_finder(double input_data [], int size){
	std::vector<double> moments(4,0.0);
	int num_elem = size;
	double n = 0.0;
	double inv_n = 0.0;
	double val = 0.0;
	double delta = 0.0;
	double A = 0.0;
	double B = 0.0;
	double mean = 0.0;
	double mom2 = 0.0;
	double mom3 = 0.0;
	double mom4 = 0.0;
	for (int i = 0; i != num_elem; ++i){
		n = i + 1.0;
		inv_n = 1.0 / n;
		val = input_data[i];
		delta = val - mean;
		A = delta * inv_n;
		mean += A;
		mom4 += A * (A * A * delta * i * (n * (n - 3.) + 3.) + 6. * A * mom2 - 4. * mom3);
		B = val - mean;
    mom3 += A * (B * delta * (n - 2.) - 3. * mom2);
    mom2 += delta * B;
	}
	moments.at(0) = mean;
	moments.at(1) = mom2;
	moments.at(2) = mom3;
	moments.at(3) = mom4;
	return moments;
}
std::vector<double> initial_moment_finder(const std::vector<double> & input_data, int size){
	std::vector<double> moments(4,0.0);
	int num_elem = size;
	double n = 0.0;
	double inv_n = 0.0;
	double val = 0.0;
	double delta = 0.0;
	double A = 0.0;
	double B = 0.0;
	double mean = 0.0;
	double mom2 = 0.0;
	double mom3 = 0.0;
	double mom4 = 0.0;
	for (int i = 0; i != num_elem; ++i){
		n = i + 1.0;
		inv_n = 1.0 / n;
		val = input_data.at(i);
		delta = val - mean;
		A = delta * inv_n;
		mean += A;
		mom4 += A * (A * A * delta * i * (n * (n - 3.) + 3.) + 6. * A * mom2 - 4. * mom3);
		B = val - mean;
    mom3 += A * (B * delta * (n - 2.) - 3. * mom2);
    mom2 += delta * B;
	}
	moments.at(0) = mean;
	moments.at(1) = mom2;
	moments.at(2) = mom3;
	moments.at(3) = mom4;
	return moments;
}
std::vector<double> moment_updater(double old_val, double new_val, int size, std::vector<double> initial_moments){
	double nm1 = size - 1.0;
	double nm2 = size - 2.0;
	double old_val_min_mu = old_val - initial_moments[0];
	double mom2 = initial_moments[1];
	double mom3 = initial_moments[2];
	double mom4 = initial_moments[3];
	double mean = 0.0;
	double nm1_sqr = nm1 * nm1;
	double old_val_min_mu_sqr = old_val_min_mu * old_val_min_mu;
	// calculate 2nd moment with old value removed
	mom2 -= size*old_val_min_mu_sqr / nm1;
	// calculate 3rd moment with old value removed
	mom3 -= old_val_min_mu * (nm2 * size * old_val_min_mu_sqr - 3 * mom2 * nm1) / (nm1_sqr);
	// calculate 4th moment with old value removed
	mom4 -= old_val_min_mu * (6 * mom2 * nm1 * old_val_min_mu - 4 * mom3 * nm1_sqr + size * (size * size - 3 * size + 3)*old_val_min_mu_sqr * old_val_min_mu) / (nm1_sqr * nm1);
	// calculate the mean with old value removed
	mean = (size * initial_moments[0] - old_val) / nm1;
	// update the moments with the new value
	double delta = new_val - mean;
	double inv_n = 1.0 / size;
	double A = delta * inv_n;
	mean += A;
	mom4 += A * (A * A * delta * nm1 * (size * (size - 3.) + 3.) + 6. * A * mom2 - 4. * mom3);
	double B = new_val - mean;
	mom3 += A * (B * delta * nm2 - 3. * mom2);
	mom2 += delta * B;
	std::vector<double> new_moments(4,0.0);
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
std::map<std::string, std::vector<double>> vs_file_reader(const char* filename){
	//open file up for reading
	std::string file_path(filename);
	std::ifstream in_s(file_path);
	std::string line;
	std::vector<std::string> line_split;
	std::map<std::string, std::vector<double>> chr_data;
	// read the header in and parse the header
	if (in_s){
		getline(in_s, line);
		std::unordered_map<std::string, int> header_indexes = header_parser(line, '\t'); // O(1) access
		// read each line, split it, check if chromosome name already exists in a set. If so, add to current vector. If new chromosome name, push current vector into map, clear the vector, add the new name to the set,
		// and add new value
		int chr_index = header_indexes.at("Chr");
		int data_index = header_indexes.at("Data");
		std::unordered_set<std::string> chr_names;
		std::vector<double> cur_chr_vals;
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
					chr_data.insert(std::pair<std::string, std::vector<double>> (old_chr_name, cur_chr_vals));
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
		chr_data.insert(std::pair<std::string, std::vector<double>> (old_chr_name, cur_chr_vals));
		in_s.close();
	}
	else {
		std::cout << "Could not open: " << file_path << std::endl;
	}
	return chr_data;
}
std::vector<double> discriminant_values_generator(const std::vector<double> & input_data, int window_size){
	int slide_iterations = input_data.size() - window_size;
	std::vector<double> initial_datapoints(window_size, 0.0);
	for (int i = 0; i != window_size; ++i){
		initial_datapoints.at(i) = input_data.at(i);
	}
	std::vector<double> moments = initial_moment_finder(initial_datapoints, window_size);


	// calculate the first discrimanant value and store it into a vector of discrimanant values
	std::vector<double> discriminants(slide_iterations + 1, 0.0);
	discriminants.at(0) = discriminant_function(200.0, 0.05, moments.at(1), moments.at(3), window_size);


	// implement the sliding window
	int old_value_index = 0;
	int new_value_index = window_size;
	//std::vector<double> new_moments;
	for (int i = 0; i != slide_iterations; ++i){
		// first calculate new moments
		moments = moment_updater(input_data[old_value_index], input_data[new_value_index], window_size, moments);
		// overwrite the old moments with the new ones
		//moments = new_moments;
		// calculate the new discrimanant and add to vector
		discriminants.at(i + 1) = discriminant_function(200.0, 0.05, moments.at(1), moments.at(3), window_size);
		++old_value_index;
		++new_value_index;
	}
	return discriminants;
}

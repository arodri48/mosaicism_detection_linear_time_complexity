import numpy as np
import sandia_stats
from collections import Counter
from math import floor

class Child:

    def __init__(self, child_name, father_name, mother_name):
        self.name = child_name
        self.father_name = father_name
        self.mother_name = mother_name
        self.pos_arr = np.empty(1)
        self.dad_rd_array = np.empty(1)
        self.mom_rd_array = np.empty(1)
        self.vcf_pos_start_of_mosaicism = 0
        self.index_diff_arr_start_of_mosaicism = 0
        self.vcf_pos_end_of_mosaicism = 0
        self.index_diff_arr_end_of_mosaicism = 0
        self.left_border_mosaicism_region = 0
        self.right_border_mosaicism_region = 0
        self.t_values = []
        self.naive_t_values = []
        self.forward_filter_results = []
        self.backward_filter_results = []
        self.precise_index_start_mosaicism = 0
        self.precise_index_end_mosaicism = 0
        self.precise_vcf_start_pos = 0
        self.precise_vcf_end_pos = 0

    def phasable_snp_determiner(self, chr_df):
        # first make temporary helper variables
        pos_final = []
        dad_rd_final = []
        mom_rd_final = []
        het_set = {"0/1", "1/0", "1|0", "0|1"}
        # iterate through every row
        for index, row in chr_df.iterrows():
            child_info = row[self.name].split(':', 3)
            if child_info[0] in het_set:
                child_read_depths = child_info[1].split(',')
                # print(child_read_depths)
                child_rd_first = int(child_read_depths[0])
                child_rd_second = int(child_read_depths[1])
                if (4 < child_rd_first < 75) and (4 < child_rd_second < 75):
                    # proband is a het; need to check the parents and check at least one is homozygous and if they are
                    # both homozygous, not for same allele
                    mom_line_info = row[self.mother_name].split(':', 2)
                    dad_line_info = row[self.father_name].split(':', 2)
                    mom_geno_count = Counter(mom_line_info[0])
                    dad_geno_count = Counter(dad_line_info[0])
                    if not ((3 == len(mom_geno_count) == len(dad_geno_count)) or
                            (2 == mom_geno_count['0'] == dad_geno_count['0']) or (
                                    2 == mom_geno_count['1'] == dad_geno_count['1'])):
                        if 0 == dad_geno_count['.'] == mom_geno_count['.']:
                            # save the position number and then the read depth for the child
                            pos_final.append(row['POS'])
                            # case 1: Dad is hom var and mom is hom ref
                            if 2 == dad_geno_count['1'] == mom_geno_count['0']:
                                if child_info[0][0] == '1':
                                    dad_rd_final.append(child_rd_first)
                                    mom_rd_final.append(child_rd_second)
                                else:
                                    dad_rd_final.append(child_rd_second)
                                    mom_rd_final.append(child_rd_first)
                            # case 2: mom is hom var and dad is hom ref
                            elif 2 == dad_geno_count['0'] == mom_geno_count['1']:
                                if child_info[0][0] == '1':
                                    dad_rd_final.append(child_rd_second)
                                    mom_rd_final.append(child_rd_first)
                                else:
                                    dad_rd_final.append(child_rd_first)
                                    mom_rd_final.append(child_rd_second)
                            # case 3: Dad is a het
                            elif len(dad_geno_count) == 3:
                                # if mom is hom ref
                                if mom_geno_count['0'] == 2:
                                    if child_info[0][0] == '0':
                                        dad_rd_final.append(child_rd_second)
                                        mom_rd_final.append(child_rd_first)
                                    else:
                                        dad_rd_final.append(child_rd_first)
                                        mom_rd_final.append(child_rd_second)
                                # if mom is hom var
                                else:
                                    if child_info[0][0] == '1':
                                        dad_rd_final.append(child_rd_second)
                                        mom_rd_final.append(child_rd_first)
                                    else:
                                        dad_rd_final.append(child_rd_first)
                                        mom_rd_final.append(child_rd_second)
                            # case 4: Mom is a het
                            else:
                                # if dad is hom ref
                                if dad_geno_count['0'] == 2:
                                    if child_info[0][0] == '0':
                                        dad_rd_final.append(child_rd_first)
                                        mom_rd_final.append(child_rd_second)
                                    else:
                                        dad_rd_final.append(child_rd_second)
                                        mom_rd_final.append(child_rd_first)
                                # if dad is hom var
                                else:
                                    if child_info[0][0] == '1':
                                        dad_rd_final.append(child_rd_first)
                                        mom_rd_final.append(child_rd_second)
                                    else:
                                        dad_rd_final.append(child_rd_second)
                                        mom_rd_final.append(child_rd_first)
        self.pos_arr = pos_final
        self.mom_rd_array = np.array(mom_rd_final)
        self.dad_rd_array = np.array(dad_rd_final)

    def t_test_snps(self, samp_size=10000, t_thres=25):

        if samp_size > self.mom_rd_array.size:
            print("Sample size is larger than total number of data points")
        else:
            # Step 1: Allocate an array that will store the t
            t_values = []
            # t_values = np.empty(self.child_ref_rd.size - samp_size + 1)
            # Step 2: Generate difference array between dad and mom rd
            diff_arr = (self.dad_rd_array - self.mom_rd_array).tolist()
            # Step 3: Calculate initial moments
            moments = sandia_stats.m1_m2_moment_generator(diff_arr[0:samp_size])
            # Step 4: Calculate t-statistic for first window
            t_values.append(abs(moments[0] / (moments[1]**0.5/samp_size)))
            # Step 5: Calculate t-values for rest of positions
            counter2 = samp_size
            mom_update_func = sandia_stats.m1_m2_moment_updater
            for i in range(self.mom_rd_array.size - samp_size):
                moments = mom_update_func(moments, diff_arr[i], diff_arr[counter2], samp_size)
                counter2 += 1
                t_values.append(abs(moments[0] / (moments[1]**0.5/samp_size)))
            # Step 6: Find the t-critical value and determine if there are any samples that exceed it
            index_of_mosaicism = next((i for i, elem in enumerate(t_values) if elem > t_thres), -1)
            # if index_of_mosaicism is not equal to 0, update the start_of_mosaicism index
            self.t_values = t_values
            if index_of_mosaicism > -1:
                # obtain VCF position from pos_arr
                self.vcf_pos_start_of_mosaicism = self.pos_arr[index_of_mosaicism + samp_size -1] if index_of_mosaicism != 0 else self.pos_arr[0]
                self.index_diff_arr_start_of_mosaicism = index_of_mosaicism + samp_size - 1 if index_of_mosaicism != 0 else 0
                # a region has been found; time to find end point
                index_of_end_of_mosaicism = next(
                    (i+index_of_mosaicism+1 for i, elem in enumerate(t_values[index_of_mosaicism + 1:]) if elem < t_thres),
                    len(t_values) - 1)
                self.vcf_pos_end_of_mosaicism = self.pos_arr[index_of_end_of_mosaicism + samp_size - 1] if index_of_end_of_mosaicism != len(t_values) - 1 else self.pos_arr[-1]
                self.index_diff_arr_end_of_mosaicism = index_of_end_of_mosaicism + samp_size - 1 if index_of_end_of_mosaicism != len(t_values) - 1 else len(self.pos_arr) - 1
                output = ["Mosaicism has been detected in child ", self.name, " with approximated start and end points at VCF positions ", str(self.vcf_pos_start_of_mosaicism), " and ", str(self.vcf_pos_end_of_mosaicism), ", respectively"]
                print("".join(output))

    def naive_t_test_snps(self, samp_size):
        if samp_size > self.mom_rd_array.size:
            print("Sample size is larger than total number of data points")
        else:
            # Step 2: Generate difference array between dad and mom rd
            diff_arr = (self.dad_rd_array - self.mom_rd_array)
            # Step 3: Generate t_values naively
            t_vals = []
            for i in range(self.mom_rd_array.size - samp_size + 1):
                local_arr = diff_arr[i:samp_size + i]
                t_vals.append(abs(local_arr.mean() / (local_arr.var() / samp_size) ** 0.5))
            self.naive_t_values = t_vals

    def edge_detection(self, sample_size = 10000):
        estimated_interval_length = self.index_diff_arr_end_of_mosaicism - self.index_diff_arr_start_of_mosaicism + 1
        width_of_average = floor(estimated_interval_length/2)
        fourth_up = floor(estimated_interval_length/4)
        diff_arr = self.dad_rd_array - self.mom_rd_array
        height = diff_arr[self.index_diff_arr_start_of_mosaicism + fourth_up:self.index_diff_arr_start_of_mosaicism + fourth_up + width_of_average].mean()
        filter_width_one_side = floor(0.25 * sample_size)
        forward_filter = np.zeros(2*filter_width_one_side)
        forward_filter[filter_width_one_side:] = height
        backward_filter = np.zeros(2*filter_width_one_side)
        backward_filter[:filter_width_one_side] = height

        if self.index_diff_arr_start_of_mosaicism == 0 and self.index_diff_arr_end_of_mosaicism == len(self.pos_arr) - 1:
            # entire chromosome is mosaic
            self.precise_index_start_mosaicism = 0
            self.precise_index_end_mosaicism = len(self.pos_arr) - 1
            self.precise_vcf_start_pos = self.pos_arr[0]
            self.precise_vcf_end_pos = self.pos_arr[-1]
            # return [self.pos_arr[0], self.pos_arr[-1]]
        elif self.index_diff_arr_start_of_mosaicism == 0 and self.index_diff_arr_end_of_mosaicism != len(self.pos_arr) - 1:
            # mosaic region from start of chromosome up to somewhere in the middle
            center_index = self.index_diff_arr_end_of_mosaicism - floor(0.5 * sample_size)
            filter_difference = [abs((diff_arr[center_index - filter_width_one_side + i: center_index + filter_width_one_side + i] - backward_filter).sum(dtype=float)) for i in range(sample_size)]
            min_val = min(filter_difference)

            self.precise_index_start_mosaicism = 0
            self.precise_index_end_mosaicism = center_index + filter_difference.index(min_val)
            self.precise_vcf_start_pos = self.pos_arr[0]
            self.precise_vcf_end_pos = self.pos_arr[self.precise_index_end_mosaicism]

            # return [self.pos_arr[0], self.pos_arr[center_index + filter_difference.index(min_val)]]
        elif self.index_diff_arr_start_of_mosaicism != 0 and self.index_diff_arr_end_of_mosaicism == len(self.pos_arr) - 1:
            # mosaic region starts somewhere in the middle and goes up to the end
            center_index = self.index_diff_arr_start_of_mosaicism - floor(0.5 * sample_size)
            filter_difference = [abs((diff_arr[center_index - filter_width_one_side + i: center_index + filter_width_one_side + i] - forward_filter).sum(dtype=float)) for i in range(sample_size)]
            min_val = min(filter_difference)

            self.precise_index_start_mosaicism = center_index + filter_difference.index(min_val)
            self.precise_index_end_mosaicism = len(self.pos_arr) - 1
            self.precise_vcf_start_pos = self.pos_arr[self.precise_index_start_mosaicism]
            self.precise_vcf_end_pos = self.pos_arr[-1]

            # return [self.pos_arr[center_index + filter_difference.index(min_val)], self.pos_arr[-1]]
        else:
            # first find the left boundary
            center_start_index = self.index_diff_arr_start_of_mosaicism - floor(0.5 * sample_size)
            filter_start_difference = [abs((diff_arr[
                                      center_start_index - filter_width_one_side + i: center_start_index + filter_width_one_side + i] - forward_filter).sum(
                dtype=float)) for i in range(sample_size)]
            min_start_val = min(filter_start_difference)
            # then find the right boundary
            center_end_index = self.index_diff_arr_end_of_mosaicism - floor(0.5 * sample_size)
            filter_end_difference = [abs((diff_arr[
                                      center_end_index - filter_width_one_side + i: center_end_index + filter_width_one_side + i] - backward_filter).sum(
                dtype=float)) for i in range(sample_size)]
            min_end_val = min(filter_end_difference)

            # save the indices and start and end positions
            self.precise_index_start_mosaicism = center_start_index + filter_start_difference.index(min_start_val)
            self.precise_index_end_mosaicism = center_end_index + filter_end_difference.index(min_end_val)
            self.precise_index_start_mosaicism = self.pos_arr[self.precise_index_start_mosaicism]
            self.precise_index_end_mosaicism = self.pos_arr[self.precise_index_end_mosaicism]
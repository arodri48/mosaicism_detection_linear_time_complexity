from collections import Counter
import numpy as np
import sandia_stats
import helper_functions
from math import floor


def sandia_t_test_snps(vcf_pos, maternal_rd, paternal_rd, samp_size=10000, t_thres=25):
    mom_minus_samp_size = maternal_rd.size - samp_size
    if 0 > mom_minus_samp_size:
        print("Sample size is larger than total number of data points")
        return None
    else:
        # Step 1: Allocate an array that will store the t
        t_values = []
        # Step 2: Generate difference array between mom and dad
        diff_arr = (maternal_rd - paternal_rd).tolist()
        # Step 3: Calculate initial moments
        moments = sandia_stats.m1_m2_moment_generator(diff_arr[0:samp_size])
        # Step 4: Calculate t-statistic for first window
        t_values.append(abs(moments[0] / (moments[1] ** 0.5 / samp_size)))
        # Step 5: Calculate t-values for rest of positions
        counter2 = samp_size
        mom_update_func = sandia_stats.m1_m2_moment_updater
        for i in range(mom_minus_samp_size):
            moments = mom_update_func(moments, diff_arr[i], diff_arr[counter2], samp_size)
            counter2 += 1
            t_values.append(abs(moments[0] / (moments[1] ** 0.5 / samp_size)))
        # Step 6: See if any t-values exceed the t-value threshold
        index_of_mosaicism = next((i for i, elem in enumerate(t_values) if elem > t_thres), -1)
        if index_of_mosaicism != -1:
            # there is a t_value that exceeds the thresold, save the vcf position of the start
            vcf_pos_start_of_mosaicism = vcf_pos[index_of_mosaicism + samp_size - 1] if index_of_mosaicism != 0 else \
                vcf_pos[0]
            # figure out the end point of mosaicism
            index_of_end_of_mosaicism = next(
                (i + index_of_mosaicism + 1 for i, elem in enumerate(t_values[index_of_mosaicism + 1:]) if
                 elem < t_thres), len(t_values) - 1)
            vcf_pos_end_of_mosaicism = vcf_pos[
                index_of_end_of_mosaicism + samp_size - 1] if index_of_end_of_mosaicism != len(t_values) - 1 else \
                vcf_pos[-1]
            return [vcf_pos_start_of_mosaicism, vcf_pos_end_of_mosaicism]
        else:
            return None


def phasable_snp_determiner(chr_df, proband_name, father_name, mother_name):
    # first make temporary helper variables
    pos_final = []
    dad_rd_final = []
    mom_rd_final = []
    het_set = {"0/1", "1/0", "1|0", "0|1"}
    # iterate through every row
    for index, row in chr_df.iterrows():
        child_info = row[proband_name].split(':', 3)
        if child_info[0] in het_set:
            child_read_depths = child_info[1].split(',')
            child_rd_first = int(child_read_depths[0])
            child_rd_second = int(child_read_depths[1])
            if (4 < child_rd_first < 75) and (4 < child_rd_second < 75):
                # proband is a het; need to check the parents and check at least one is homozygous and if they are
                # both homozygous, not for same allele
                mom_line_info = row[mother_name].split(':', 2)
                dad_line_info = row[father_name].split(':', 2)
                mom_geno_count = Counter(mom_line_info[0])
                dad_geno_count = Counter(dad_line_info[0])
                if not (dad_geno_count['.'] or mom_geno_count['.']):
                    is_mom_hom_ref = 2 == mom_geno_count['0']
                    is_dad_hom_ref = 2 == dad_geno_count['0']
                    is_mom_hom_var = 2 == mom_geno_count['1']
                    is_dad_hom_var = 2 == dad_geno_count['1']
                    is_dad_het = 3 == len(dad_geno_count)
                    is_mom_het = 3 == len(mom_geno_count)
                    if not ((is_dad_het and is_mom_het) or
                            (is_mom_hom_ref and is_dad_hom_ref) or (
                                    is_mom_hom_var and is_dad_hom_var)):
                        # save the position number and then the read depth for the child
                        pos_final.append(row['POS'])
                        # case 1: Dad is hom var and mom is hom ref
                        if is_dad_hom_var and is_mom_hom_ref:
                            if child_info[0][0] == '1':
                                dad_rd_final.append(child_rd_first)
                                mom_rd_final.append(child_rd_second)
                            else:
                                dad_rd_final.append(child_rd_second)
                                mom_rd_final.append(child_rd_first)
                        # case 2: mom is hom var and dad is hom ref
                        elif is_mom_hom_var and is_dad_hom_ref:
                            if child_info[0][0] == '1':
                                dad_rd_final.append(child_rd_second)
                                mom_rd_final.append(child_rd_first)
                            else:
                                dad_rd_final.append(child_rd_first)
                                mom_rd_final.append(child_rd_second)
                        # case 3: Dad is a het
                        elif is_dad_het:
                            # if mom is hom ref
                            if is_mom_hom_ref:
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
                            if is_dad_hom_ref:
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
    return pos_final, np.array(mom_rd_final), np.array(dad_rd_final)


def edge_detection(sample_size, estimated_start_index, estimated_end_index, paternal_rd_array, maternal_rd_array):
    estimated_interval_length = estimated_end_index - estimated_start_index
    width_of_average = floor(estimated_interval_length / 2)
    fourth_up = floor(estimated_interval_length / 4)
    diff_arr = paternal_rd_array - maternal_rd_array
    height = diff_arr[estimated_start_index + fourth_up:estimated_start_index + fourth_up + width_of_average].mean()
    filter_width_one_side = floor(0.25 * sample_size)
    forward_filter = np.zeros(2 * filter_width_one_side)
    forward_filter[filter_width_one_side:] = height
    backward_filter = np.zeros(2 * filter_width_one_side)
    backward_filter[:filter_width_one_side] = height

    final_index = paternal_rd_array.size - 1
    is_mosaicism_to_the_end = estimated_end_index == final_index

    if not estimated_start_index and is_mosaicism_to_the_end:
        # estimated chromosome is mosaic
        return [0, final_index]
    elif not (estimated_start_index or is_mosaicism_to_the_end):
        center_index = estimated_end_index - floor(0.5 * sample_size)
        filter_difference = [abs((diff_arr[
                                  center_index - filter_width_one_side + i: center_index + filter_width_one_side + i] - backward_filter).sum(
            dtype=float)) for i in range(sample_size)]
        min_val = min(filter_difference)
        return [0, center_index + filter_difference.index(min_val)]
    elif estimated_start_index and is_mosaicism_to_the_end:
        center_index = estimated_start_index - floor(0.5 * sample_size)
        filter_difference = [abs((diff_arr[
                                  center_index - filter_width_one_side + i: center_index + filter_width_one_side + i] - forward_filter).sum(
            dtype=float)) for i in range(sample_size)]
        min_val = min(filter_difference)
        return [center_index + filter_difference.index(min_val), final_index]
    else:
        center_start_index = estimated_start_index - floor(0.5 * sample_size)
        filter_start_difference = [abs((diff_arr[
                                        center_start_index - filter_width_one_side + i: center_start_index + filter_width_one_side + i] - forward_filter).sum(
            dtype=float)) for i in range(sample_size)]
        min_start_val = min(filter_start_difference)

        center_end_index = estimated_end_index - floor(0.5 * sample_size)
        filter_end_difference = [abs((diff_arr[
                                      center_end_index - filter_width_one_side + i: center_end_index + filter_width_one_side + i] - backward_filter).sum(
            dtype=float)) for i in range(sample_size)]
        min_end_val = min(filter_end_difference)

        return [center_start_index + filter_start_difference.index(min_start_val),
                center_end_index + filter_end_difference.index(min_end_val)]


def runner(child, chr_name, sample_size, t_threshold, SNP_df):
    # step 1: filter the SNP df by chromosome name
    chr_snp_df = helper_functions.chromosome_filter(SNP_df, chr_name)
    # step 2: do phasing and return results
    vcf_pos, maternal_rd, paternal_rd = phasable_snp_determiner(chr_snp_df, child.name, child.father_name,
                                                                child.mother_name)
    mosaicism_initial_survey_results = sandia_t_test_snps(vcf_pos, maternal_rd, paternal_rd, samp_size=sample_size, t_thres=t_threshold)
    # TODO: change output of sandia t test snps to make sure program works properly (have it return approximate indices instead of vcf_pos)
    # TODO: then have edge detection return the VCF pos
    return edge_detection(sample_size, mosaicism_initial_survey_results[0], mosaicism_initial_survey_results[1], paternal_rd, maternal_rd) if mosaicism_initial_survey_results is not None else None

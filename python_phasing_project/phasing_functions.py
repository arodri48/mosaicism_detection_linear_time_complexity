from collections import Counter
import numpy as np
import sandia_stats
from math import floor
from scipy import stats


def sandia_t_test_snps(maternal_rd, paternal_rd, samp_size=10000, t_thres=25):
    mom_minus_samp_size = maternal_rd.size - samp_size
    # Check if sample size is greater read depth total
    if 0 > mom_minus_samp_size:
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
            index_start_of_mosaicism_vcf = index_of_mosaicism + samp_size - 1 if index_of_mosaicism != 0 else 0
            # figure out the end point of mosaicism
            index_of_end_of_mosaicism = next(
                (i + index_of_mosaicism + 1 for i, elem in enumerate(t_values[index_of_mosaicism + 1:]) if
                 elem < t_thres), len(t_values) - 1)
            index_end_of_mosaicism_vcf = index_of_end_of_mosaicism + samp_size - 1
            return [index_start_of_mosaicism_vcf, index_end_of_mosaicism_vcf]
        else:
            return None


def phasable_snp_determiner(chr_df, proband_name, father_name, mother_name):
    # first make temporary helper variables
    pos_final = []
    dad_rd_final = []
    mom_rd_final = []
    father_total_read_depth = []
    mother_total_read_depth = []
    het_set = {"0/1", "1/0", "1|0", "0|1"}

    for row in chr_df.itertuples(index=False):
        child_info = getattr(row, proband_name).split(':', 2)
        if child_info[0] in het_set:
            mom_line_info = getattr(row, mother_name).split(':', 2)
            dad_line_info = getattr(row, father_name).split(':', 2)
            mom_geno_count = Counter(mom_line_info[0])
            dad_geno_count = Counter(dad_line_info[0])
            # quickly check if both parents are het, if so, not phasable
            is_dad_het = 3 == len(dad_geno_count)
            is_mom_het = 3 == len(mom_geno_count)
            if not (is_dad_het and is_mom_het):
                if not (dad_geno_count['.'] or mom_geno_count['.']):
                    is_mom_hom_ref = 2 == mom_geno_count['0']
                    is_dad_hom_ref = 2 == dad_geno_count['0']
                    is_mom_hom_var = 2 == mom_geno_count['1']
                    is_dad_hom_var = 2 == dad_geno_count['1']
                    if not ((is_mom_hom_ref and is_dad_hom_ref) or (is_mom_hom_var and is_dad_hom_var)):
                        child_read_depths = child_info[1].split(',', 2)
                        child_rd_first = int(child_read_depths[0])
                        child_rd_second = int(child_read_depths[1])
                        if (4 < child_rd_first < 75) and (4 < child_rd_second < 75):
                            dad_read_depths = dad_line_info[1].split(',', 2)
                            dad_rd_first = int(dad_read_depths[0])
                            dad_rd_second = int(dad_read_depths[1])
                            mom_read_depths = mom_line_info[1].split(',', 2)
                            mom_rd_first = int(mom_read_depths[0])
                            mom_rd_second = int(mom_read_depths[1])
                            # save the position number and then the read depth for the child
                            pos_final.append(getattr(row, 'POS'))
                            father_total_read_depth.append(dad_rd_first + dad_rd_second)
                            mother_total_read_depth.append(mom_rd_first + mom_rd_second)
                            # case 1: Dad is a het
                            if is_dad_het:
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
                            # case 2: Mom is a het
                            elif is_mom_het:
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
                            # case 3: dad is hom var and mom is hom ref
                            elif is_dad_hom_var and is_mom_hom_ref:
                                if child_info[0][0] == '1':
                                    dad_rd_final.append(child_rd_first)
                                    mom_rd_final.append(child_rd_second)
                                else:
                                    dad_rd_final.append(child_rd_second)
                                    mom_rd_final.append(child_rd_first)
                            # case 4: mom is hom var and dad is hom ref
                            else:
                                if child_info[0][0] == '1':
                                    dad_rd_final.append(child_rd_second)
                                    mom_rd_final.append(child_rd_first)
                                else:
                                    dad_rd_final.append(child_rd_first)
                                    mom_rd_final.append(child_rd_second)

    return pos_final, np.array(mom_rd_final), np.array(dad_rd_final), np.array(father_total_read_depth), np.array(mother_total_read_depth)


def edge_detection(sample_size, estimated_start_index, estimated_end_index, paternal_rd_array, maternal_rd_array,
                   vcf_pos, edge_detection_width):
    estimated_interval_length = estimated_end_index - estimated_start_index
    width_of_average = floor(estimated_interval_length / 2)
    fourth_up = floor(estimated_interval_length / 4)
    diff_arr = paternal_rd_array - maternal_rd_array
    height = diff_arr[estimated_start_index + fourth_up:estimated_start_index + fourth_up + width_of_average].mean()
    filter_width_one_side = floor(edge_detection_width / 2 * sample_size)
    forward_filter = np.zeros(2 * filter_width_one_side)
    forward_filter[filter_width_one_side:] = height
    backward_filter = np.zeros(2 * filter_width_one_side)
    backward_filter[:filter_width_one_side] = height

    final_index = paternal_rd_array.size - 1
    is_mosaicism_to_the_end = estimated_end_index == final_index
    # case 1: whole chromosome is mosaic
    if not estimated_start_index and is_mosaicism_to_the_end:
        return [vcf_pos[0], vcf_pos[final_index], 0, final_index, height]
    # case 2: mosaicism starts and beginning and ends in the middle
    elif not (estimated_start_index or is_mosaicism_to_the_end):
        center_index = estimated_end_index - floor(0.5 * sample_size)
        filter_difference = [abs((diff_arr[
                                  center_index - filter_width_one_side + i: center_index + filter_width_one_side + i] - backward_filter).sum(
            dtype=float)) for i in range(sample_size)]
        min_val = min(filter_difference)
        return [vcf_pos[0], vcf_pos[center_index + filter_difference.index(min_val)], 0,
                center_index + filter_difference.index(min_val), height]
    # case 3: mosaicism starts in middle and ends at the telomere
    elif estimated_start_index and is_mosaicism_to_the_end:
        center_index = estimated_start_index - floor(0.5 * sample_size)
        filter_difference = [abs((diff_arr[
                                  center_index - filter_width_one_side + i: center_index + filter_width_one_side + i] - forward_filter).sum(
            dtype=float)) for i in range(sample_size)]
        min_val = min(filter_difference)
        return [vcf_pos[center_index + filter_difference.index(min_val)], vcf_pos[final_index],
                center_index + filter_difference.index(min_val), final_index, height]
    # case 4: mosaicism starts and ends somewhere in teh middle
    else:
        center_start_index = estimated_start_index - floor(0.5 * sample_size)
        filter_start_difference = [abs((diff_arr[
                                        center_start_index - filter_width_one_side + i: center_start_index + filter_width_one_side + i] - forward_filter).sum(
            dtype=float)) for i in range(sample_size)]
        min_start_val = min(filter_start_difference)

        center_end_index = estimated_end_index - floor(0.5 * sample_size)
        filter_end_difference = [abs((diff_arr[
                                      center_end_index - filter_width_one_side + i:center_end_index + filter_width_one_side + i] - backward_filter).sum(
            dtype=float)) for i in range(sample_size)]
        min_end_val = min(filter_end_difference)

        return [vcf_pos[center_start_index + filter_start_difference.index(min_start_val)],
                vcf_pos[center_end_index + filter_end_difference.index(min_end_val)],
                center_start_index + filter_start_difference.index(min_start_val),
                center_end_index + filter_end_difference.index(min_end_val), height]


def runner(child_name, father_name, mother_name, sample_size, t_threshold, chr_snp_df, edge_detection_width):
    # step 1: do phasing
    vcf_pos, maternal_rd, paternal_rd, father_total_rd, mother_total_rd = phasable_snp_determiner(chr_snp_df, child_name, father_name, mother_name)
    # step 2: determine if mosaicism is present
    mosaicism_initial_survey_results = sandia_t_test_snps(maternal_rd, paternal_rd, samp_size=sample_size,
                                                          t_thres=t_threshold)
    # step 3: if it is, accurately determine edges
    if mosaicism_initial_survey_results is None:
        return None
    else:
        edge_detection_results = edge_detection(sample_size, mosaicism_initial_survey_results[0],
                                                mosaicism_initial_survey_results[1],
                                                paternal_rd, maternal_rd, vcf_pos,
                                                edge_detection_width)
        y_mosaicism_detector_results = y_mosaicism_detector(edge_detection_results, maternal_rd, paternal_rd)
        if y_mosaicism_detector_results[0] == 1:
            return edge_detection_results + y_mosaicism_detector_results
        else:
            quantification_results = mosaicism_quantifier(maternal_rd, paternal_rd, edge_detection_results[2], edge_detection_results[3], father_total_rd, mother_total_rd)
            return edge_detection_results + y_mosaicism_detector_results + quantification_results + [maternal_rd, paternal_rd]

def mosaicism_quantifier(mat_rd, pat_rd, start_region_index, end_region_index, father_total_rd, mother_total_rd):
    # Step 1: Extract mosaic region
    mosaic_region_dad = pat_rd[start_region_index:end_region_index]
    mosaic_region_mom = mat_rd[start_region_index:end_region_index]
    # Step 2: Calculate allele depth ratio
    dad_rd_sum = mosaic_region_dad.sum()
    mom_rd_sum = mosaic_region_mom.sum()
    mat_allele_depth_ratio = mom_rd_sum / (mom_rd_sum + dad_rd_sum)
    # Step 4: Obtain fraction of mosaicism
    mosaic_percentage = mosaicism_quantification(mat_allele_depth_ratio)
    return mosaic_percentage

def mosaicism_quantification(allele_depth_ratio):
    # type number 1: monosomy-disomy
    # type number 2: trisomy-disomy
    # type number 3: UPD-disomy

    if 0.34 < allele_depth_ratio < 0.5:
        return [(2 * allele_depth_ratio - 1) / (allele_depth_ratio - 1), (1 / allele_depth_ratio) - 2, 1 - 2 * allele_depth_ratio]
    elif allele_depth_ratio <= 0.34:
        return [(2 * allele_depth_ratio - 1) / (allele_depth_ratio - 1), 0.95, 1 - 2 * allele_depth_ratio]
    elif 0.5 < allele_depth_ratio < 0.66:
        return [2 - (1 / allele_depth_ratio), (2 * allele_depth_ratio - 1) / (1 - allele_depth_ratio), 2 * allele_depth_ratio - 1]
    else:
        return [2 - (1 / allele_depth_ratio), 0.95, 2 * allele_depth_ratio - 1]


def y_mosaicism_detector(edge_detection_results, mat_rd, pat_rd):
    # Step 1: save the start and end indices
    start = edge_detection_results[2]
    end = edge_detection_results[3]
    mosaic_range = end - start
    # Step 2: Obtain middle 50% of read depth difference
    diff_arr = pat_rd - mat_rd
    mid_mosaic_region = diff_arr[start + floor(0.25 * mosaic_range): end - floor(0.25 * mosaic_range)]
    # Step 3: Do a t-test to see if slope is non-zero
    x_axis = np.array([i for i in range(start + floor(0.25 * mosaic_range), end - floor(0.25 * mosaic_range))])
    res = stats.linregress(x_axis, mid_mosaic_region)
    # Step 4: Check result
    if abs(res.slope) < 0.01:
        return [0]
    else:
        return [1, res.slope, floor(-1 * res.intercept / res.slope)]


def t_test_runner(child_name, father_name, mother_name, sample_size, chr_snp_df):
    # step 1: do phasing
    vcf_pos, maternal_rd, paternal_rd = phasable_snp_determiner(chr_snp_df, child_name, father_name, mother_name)
    t_test_vals = no_classifier_t_test(maternal_rd, paternal_rd, samp_size=sample_size)
    return t_test_vals


def no_classifier_t_test(maternal_rd, paternal_rd, samp_size=10000):
    mom_minus_samp_size = maternal_rd.size - samp_size
    # Check if sample size is greater read depth total
    if 0 > mom_minus_samp_size:
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
        return t_values

def mosaicism_classification_function(mosaicism_outcome, vcf_file, dad_name, mom_name):
    # Step 4a: First check if any are not None
    any_mosaicism = any(elem is not None for elem in mosaicism_outcome)
    if not any_mosaicism:
        return None
    else:
        # There is at least one detected area of mosaicism
        # Step 4b: If so, generate a list of 22 zeros
        classification_results = [0 for elem in range(22)]
        # Step 4c: Generate a dictionary of VCF lines by chromosome
        vcf_file_by_chr = dict(tuple(vcf_file.groupby(['#CHROM'])))
        # Step 4d: Obtain average read depth for each parent for each chromosome
        read_depth_averages = [avg_read_depth_finder(vcf_file_by_chr[chr_name], mom_name, dad_name) for chr_name in range(1, 23)]
        # Step 4e: Obtain overall average read depth per parent
        dad_avg_read_depths = [elem[0] for elem in read_depth_averages]
        mom_avg_read_depths = [elem[1] for elem in read_depth_averages]
        dad_grand_avg = sum(dad_avg_read_depths)/len(dad_avg_read_depths)
        mom_grand_avg = sum(mom_avg_read_depths)/len(mom_avg_read_depths)
        # Step 4f: For each elem in mosaicism_outcome, if not None, determine most likely type of mosaicism
        for i, elem in enumerate(mosaicism_outcome):
            if elem is not None:
                # TODO: classification_results[i] = classification_determiner(heh)
        return classification_results

def classification_determiner(paternal_rd, maternal_rd, start_index, end_index, paternal_total_avg_rd, maternal_total_avg_rd):
    # I need to know the child's read depths in the mosaic region, the dad's total avg read depth, the mom's total avg
    # read depth
    dad_rd_individual_chr_avg = paternal_total_avg_rd / 2
    mom_rd_individual_chr_avg = maternal_total_avg_rd / 2
    child_paternal_rd_avg =

    # Obtain results for disomy-UDP

    # Obtain results for disomy-trisomy

    # Obtain results for disomy-monosomy

def avg_read_depth_finder(vcf_lines, mom_name, dad_name):
    dad_mean = 0.0
    mom_mean = 0.0
    n = 1
    for row in vcf_lines.itertuples(index=False):
        mom_line_info = getattr(row, mom_name).split(':', 2)
        dad_line_info = getattr(row, dad_name).split(':', 2)
        mom_geno_count = Counter(mom_line_info[0])
        dad_geno_count = Counter(dad_line_info[0])
        if not (dad_geno_count['.'] or mom_geno_count['.']):
            dad_read_depths = dad_line_info[1].split(',', 2)
            dad_rd_first = int(dad_read_depths[0])
            dad_rd_second = int(dad_read_depths[1])
            mom_read_depths = mom_line_info[1].split(',', 2)
            mom_rd_first = int(mom_read_depths[0])
            mom_rd_second = int(mom_read_depths[1])
            # increment the means
            dad_mean += (dad_rd_first + dad_rd_second - dad_mean) / n
            mom_mean += (mom_rd_first + mom_rd_second - mom_mean) / n
            n += 1
    return dad_mean, mom_mean
from collections import Counter
import numpy as np

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
        t_values.append(abs(moments[0] / (moments[1]**0.5/samp_size)))
        # Step 5: Calculate t-values for rest of positions
        counter2 = samp_size
        mom_update_func = sandia_stats.m1_m2_moment_updater
        for i in range(mom_minus_samp_size):
            moments = mom_update_func(moments, diff_arr[i], diff_arr[counter2], samp_size)
            counter2 += 1
            t_values.append(abs(moments[0] / (moments[1]**0.5/samp_size)))
        # Step 6: See if any t-values exceed the t-value threshold
        index_of_mosaicism = next((i for i, elem in enumerate(t_values) if elem > t_thres), -1)
        if index_of_mosaicism != -1:
            # there is a t_value that exceeds the thresold, save the vcf position of the start
            vcf_pos_start_of_mosaicism = vcf_pos[index_of_mosaicism + samp_size -1] if index_of_mosaicism != 0 else vcf_pos[0]
            # figure out the end point of mosaicism
            index_of_end_of_mosaicism = next((i+index_of_mosaicism+1 for i, elem in enumerate(t_values[index_of_mosaicism + 1:]) if elem < t_thres),len(t_values) - 1)
            vcf_pos_end_of_mosaicism = vcf_pos[index_of_end_of_mosaicism + samp_size - 1] if index_of_end_of_mosaicism != len(t_values) - 1 else vcf_pos[-1]
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
            # print(child_read_depths)
            child_rd_first = int(child_read_depths[0])
            child_rd_second = int(child_read_depths[1])
            if (4 < child_rd_first < 75) and (4 < child_rd_second < 75):
                # proband is a het; need to check the parents and check at least one is homozygous and if they are
                # both homozygous, not for same allele
                mom_line_info = row[mother_name].split(':', 2)
                dad_line_info = row[father_name].split(':', 2)
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
    return pos_final, np.array(mom_rd_final), np.array(dad_rd_final)

    
def runner(child, chr_name, sample_size, t_threshold, SNP_df):
    # step 1: filter the SNP df by chromosome name
    chr_snp_df = helper_functions.chromosome_filter(SNP_df, chr_name)
    # step 2: do phasing and return results
    vcf_pos, maternal_rd, paternal_rd = phasable_snp_determiner(chr_snp_df, child.name, child.father_name, child.mother_name)
    results = sandia_t_test_snps(vcf_pos, maternal_rd, paternal_rd, samp_size=sample_size, t_thres=t_threshold)
    return results
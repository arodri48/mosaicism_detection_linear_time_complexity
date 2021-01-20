import numpy as np
import sandia_stats
from scipy.stats import t


class Child:


    def __init__(self, child_name, father_name, mother_name):
        self.name = child_name
        self.father_name = father_name
        self.mother_name = mother_name
        self.pos_arr = np.zeros(2)
        self.dad_rd_array = np.zeros(2)
        self.mom_rd_array = np.zeros(2)
        self.dad_ref_var_array = np.zeros(2)
        self.mom_ref_var_array = np.zeros(2)
        self.child_ref_rd = np.zeros(2)
        self.child_var_rd = np.zeros(2)
        self.est_start_of_mosaicism = 0

    def phasable_snp_determiner(self, chr_df):
        # first make temporary helper variables
        pos_final = []
        dad_ref_var_final = []
        mom_ref_var_final = []
        dad_rd_final = []
        mom_rd_final = []
        child_ref_rd_final = []
        child_var_rd_final = []
        het_set = {"0/1", "1/0", "1|0", "0|1"}
        hom_ref = {"0/0", "0|0"}
        hom_var = {"1/1", "1|1"}
        het_zero_one = {"0/1", "0|1"}
        empty_geno = {"./."}

        # iterate through every row
        for index, row in chr_df.iterrows():
            child_info = row[self.name].split(':')
            child_read_depths = child_info[1].split(',')
            if child_info[0] in het_set and (4 < int(child_read_depths[0]) < 75) and (4 < int(child_read_depths[1]) < 75):
                # proband is a het; need to check the parents and check at least one is homozygous and if they are
                # both homozygous, not for same allele
                mom_line_info = row[self.mother_name].split(':')
                dad_line_info = row[self.father_name].split(':')
                dad_geno = dad_line_info[0]
                mom_geno = mom_line_info[0]
                if not((mom_geno in het_set and dad_geno in het_set) or ((mom_geno in hom_ref and dad_geno in hom_ref) or (
                        mom_geno in hom_var and dad_geno in hom_var)) ):
                    if not (dad_geno in empty_geno or mom_geno in empty_geno):
                        dad_read_depths = dad_line_info[1].split(',')
                        mom_read_depths = mom_line_info[1].split(',')
                        # save the position number and then the read depth for the child
                        pos_final.append(row['POS'])
                        if child_info[0] in het_zero_one:
                            child_ref_rd_final.append(int(child_read_depths[0]))
                            child_var_rd_final.append(int(child_read_depths[1]))
                        else:
                            child_ref_rd_final.append(int(child_read_depths[1]))
                            child_var_rd_final.append(int(child_read_depths[0]))
                        # determine which parent each allele came from and save the read depth for that parent's allele
                        if dad_geno in hom_var and mom_geno in hom_ref:
                            dad_ref_var_final.append(1)
                            mom_ref_var_final.append(0)
                            # obtain and store read count
                            dad_rd_final.append(int(dad_read_depths[0]) + int(dad_read_depths[1]))
                            mom_rd_final.append(int(mom_read_depths[0]) + int(mom_read_depths[1]))
                        elif dad_geno in hom_ref and mom_geno in hom_var:
                            dad_ref_var_final.append(0)
                            mom_ref_var_final.append(1)
                            # obtain and store read count
                            dad_rd_final.append(int(dad_read_depths[0]) + int(dad_read_depths[1]))
                            mom_rd_final.append(int(mom_read_depths[0]) + int(mom_read_depths[1]))
                        elif dad_geno in het_set:
                            if mom_geno in hom_ref:
                                mom_ref_var_final.append(0)
                                dad_ref_var_final.append(1)
                                mom_rd_final.append(int(mom_read_depths[0]) + int(mom_read_depths[1]))
                                dad_geno_split = [char for char in dad_geno]
                                if dad_geno_split[0] == '1':
                                    dad_rd_final.append(int(dad_read_depths[0]))
                                else:
                                    dad_rd_final.append(int(dad_read_depths[1]))
                            else:
                                # mom_geno is hom_var
                                mom_ref_var_final.append(1)
                                dad_ref_var_final.append(0)
                                mom_rd_final.append(int(mom_read_depths[0]) + int(mom_read_depths[1]))
                                dad_geno_split = [char for char in dad_geno]
                                if dad_geno_split[0] == '0':
                                    dad_rd_final.append(int(dad_read_depths[0]))
                                else:
                                    dad_rd_final.append(int(dad_read_depths[1]))
                        else:
                            if dad_geno in hom_ref:
                                mom_ref_var_final.append(1)
                                dad_ref_var_final.append(0)
                                dad_rd_final.append(int(dad_read_depths[0]) + int(dad_read_depths[1]))
                                mom_geno_split = [char for char in mom_geno]
                                if mom_geno_split[0] == '1':
                                    mom_rd_final.append(int(mom_read_depths[0]))
                                else:
                                    mom_rd_final.append(int(mom_read_depths[1]))
                            else:
                                mom_ref_var_final.append(0)
                                dad_ref_var_final.append(1)
                                dad_rd_final.append(int(dad_read_depths[0]) + int(dad_read_depths[1]))
                                mom_geno_split = [char for char in mom_geno]
                                if mom_geno_split[0] == '0':
                                    mom_rd_final.append(int(mom_read_depths[0]))
                                else:
                                    mom_rd_final.append(int(mom_read_depths[1]))
        self.pos_arr = np.array(pos_final)
        self.mom_rd_array = np.array(mom_rd_final)
        self.dad_rd_array = np.array(dad_rd_final)
        self.dad_ref_var_array = np.array(dad_ref_var_final)
        self.mom_ref_var_array = np.array(mom_ref_var_final)
        self.child_ref_rd = np.array(child_ref_rd_final)
        self.child_var_rd = np.array(child_var_rd_final)

    def t_test_snps(self, samp_size):

        if samp_size > self.child_ref_rd.size:
            print("Sample size is larger than total number of data points")
        else:
            nm1 = samp_size - 1
            # Step 1: Allocate an array that will store the t
            t_values = np.zeros(self.child_ref_rd.size - samp_size + 1)
            # Step 2: Generate difference array between dad and mom rd
            diff_arr = self.dad_rd_array - self.mom_rd_array
            # Step 3: Calculate initial moments
            moments = sandia_stats.statistical_moment_generator(diff_arr[0, samp_size])
            # Step 4: Calculate t-statistic for first window
            sample_var = (moments[1] / nm1)
            t_values[0] = moments[0] / ((sample_var / samp_size) ** 0.5)
            # Step 5: Calculate t-values for rest of positions
            counter1 = 0
            counter2 = samp_size
            for i in range(self.child_ref_rd.size - samp_size):
                moments = sandia_stats.moment_updater(moments, diff_arr[counter1], diff_arr[counter2], samp_size)
                counter1 += 1
                counter2 += 1
                sample_var = moments[1] / nm1
                t_values[counter1] = moments[0] / ((sample_var / samp_size) ** 0.5)
            # Step 6: Find the t-critical value and determine if there are any samples that exceed it
            t_crit = t.ppf(0.95, nm1)
            t_val_abs = np.abs(t_values)
            index_of_mosaicism = np.argmax(t_val_abs > t_crit)
            # if index_of_mosaicism is not equal to 0, update the start_of_mosaicism index
            if index_of_mosaicism > 0:
                self.est_start_of_mosaicism = index_of_mosaicism + samp_size
                print("Mosaicism has been detected in child " + self.name)

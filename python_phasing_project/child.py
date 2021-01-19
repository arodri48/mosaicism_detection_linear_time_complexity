import numpy as np
import sandia_stats


class Child:
    het_set = {"0/1", "1/0", "1|0", "0|1"}
    hom_ref = {"0/0", "0|0"}
    hom_var = {"1/1", "1|1"}
    het_zero_one = {"0/1", "0|1"}
    het_one_zero = {"1/0", "1|0"}
    empty_geno = {"./."}

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

    def phasable_snp_determiner(self, chr_df):
        # first make temporary helper variables
        mom_line_info = []
        dad_line_info = []
        child_info = []
        dad_geno = ""
        mom_geno = ""
        pos_final = []
        dad_ref_var_final = []
        mom_ref_var_final = []
        dad_rd_final = []
        mom_rd_final = []
        child_ref_rd_final = []
        child_var_rd_final = []

        # iterate through every row
        for index, row in chr_df.iterrows():
            child_info = row[self.name].split(':')
            if child_info[0] in self.het_set:
                # proband is a het; need to check the parents and check at least one is homozygous and if they are
                # both homozygous, not for same allele
                mom_line_info = row[self.mother_name].split(':')
                dad_line_info = row[self.father_name].split(':')
                dad_geno = dad_line_info[0]
                mom_geno = mom_line_info[0]
                if (not ((mom_geno in self.het_set) and (dad_geno in self.het_set)) and not (
                        (mom_geno in self.hom_ref and dad_geno in self.hom_ref) or (
                        mom_geno in self.hom_var and dad_geno in self.hom_var))):
                    if not (dad_geno in self.empty_geno or mom_geno in self.empty_geno):
                        dad_read_depths = dad_line_info[1].split(',')
                        mom_read_depths = mom_line_info[1].split(',')
                        child_read_depths = child_info[1].split(',')
                        # save the position number and then the read depth for the child
                        pos_final.append(row['POS'])
                        if child_info[0] in self.het_zero_one:
                            child_ref_rd_final.append(int(child_read_depths[0]))
                            child_var_rd_final.append(int(child_read_depths[1]))
                        else:
                            child_ref_rd_final.append(int(child_read_depths[1]))
                            child_var_rd_final.append(int(child_read_depths[0]))
                        # determine which parent each allele came from and save the read depth for that parent's allele
                        if dad_geno in self.hom_var and mom_geno in self.hom_ref:
                            dad_ref_var_final.append(1)
                            mom_ref_var_final.append(0)
                            # obtain and store read count
                            dad_rd_final.append(int(dad_read_depths[0]) + int(dad_read_depths[1]))
                            mom_rd_final.append(int(mom_read_depths[0]) + int(mom_read_depths[1]))
                        elif dad_geno in self.hom_ref and mom_geno in self.hom_var:
                            dad_ref_var_final.append(0)
                            mom_ref_var_final.append(1)
                            # obtain and store read count
                            dad_rd_final.append(int(dad_read_depths[0]) + int(dad_read_depths[1]))
                            mom_rd_final.append(int(mom_read_depths[0]) + int(mom_read_depths[1]))
                        elif dad_geno in self.het_set:
                            if mom_geno in self.hom_ref:
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
                            if dad_geno in self.hom_ref:
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
            # Step 1: Allocate an array that will store the t
            t_values = np.zeros(self.child_ref_rd.size - samp_size + 1)
            # Step 2: Generate difference array between dad and mom rd
            diff_arr = self.dad_rd_array - self.mom_rd_array
            # Step 3: Calculate initial moments
            moments = sandia_stats.statistical_moment_generator(diff_arr[0, samp_size])
            # Step 4: Calculate t-statistic for first window s

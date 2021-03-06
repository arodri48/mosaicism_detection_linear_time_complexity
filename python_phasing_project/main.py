import helper_functions
import phasing_functions


def main(config_file_path):
    # Step 0: Read config file
    with open(config_file_path, 'r') as f:
        config_elem = {}
        for line in f:
            line_split = line.strip().split('\t')
            config_elem[line_split[0]] = line_split[1]
    # Step 1: Obtain names from ped file
    names = helper_functions.ped_file_reader(config_elem["PED_FILE"])

    # Step 2: Read in the VCF file, obtain only SNPs, then group by chromosome name
    vcf_file = helper_functions.read_VCF(config_elem["VCF_FILE"], names)
    SNP_df = helper_functions.SNP_filter(vcf_file)
    snp_chr_dict = dict(tuple(SNP_df.groupby(["#CHROM"])))

    # Step 3: Generate results
    mosaicism_outcome = [phasing_functions.runner(config_elem["PROBAND_NAME"], names[0], names[1],
                                                  int(config_elem["SAMPLE_SIZE"]), float(config_elem["T_THRES"]),
                                                  snp_chr_dict[chr_name], float(config_elem["EDGE_DETECTION_WIDTH"])) for chr_name in
                         range(1, 23)]
    print(mosaicism_outcome)
    # Step 4: If any mosaicism is detected, classify it
    # Step 4a: First check if any are not None
    # Step 4b: If so, generate a list of 22 zeros
    # Step 4c: Generate a dictionary of VCF lines by chromosome
    # Step 4d: Obtain average read depth for each parent for each chromosome
    # Step 4e: Obtain overall average read depth per parent
    # Step 4f: For each elem in mosaicism_outcome, if not None, determine most likely type of mosaicism

    mosaicism_classification = phasing_functions.mosaicism_classification_function(mosaicism_outcome, vcf_file, names[0], names[1])
    print(mosaicism_classification)

    # Step 4: Write results to file
    with open(config_elem["OUTPUT_FILE"], 'w') as output_file:
        output_file.write(" ".join(["Mosaicism results for", config_elem["PROBAND_NAME"]]))
        output_file.write("\n")
        output_file.write("\t".join(["chr_number", "VCF_start", "VCF_end", "Y-Shaped_Mosaicism_present", "Slope", "Start_of_Y-Shaped_Mosaicism", "Most_Likely_Mosaic_Model", "Fraction_Under_Mono_Disomy(1)", "Fraction_Under_Tri_Disomy(2)", "Fraction_Under_UPD_Disomy(3)"]))
        output_file.write("\n")
        for i, elem in enumerate(mosaicism_outcome, start=1):
            if elem is not None:
                #print(elem)
                if elem[5] == 0:
                    output_file.write("\t".join([str(i), str(elem[0]), str(elem[1]), str(elem[5]), "N/A", "N/A", str(mosaicism_classification[i-1]), str(elem[6]), str(elem[7]), str(elem[8])]))
                else:
                    # y shaped mosaicism is present
                    output_file.write("\t".join([str(i), str(elem[0]), str(elem[1]), str(elem[5]), str(elem[6]), str(elem[7]), "N/A", "N/A", "N/A", "N/A"]))
                output_file.write("\n")



main('phasing_config_file.txt')

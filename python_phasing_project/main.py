import helper_functions
import phasing_functions
import time


def main(config_file_path):
    # Step 0: Read config file
    with open(config_file_path, 'r') as f:
        config_elem = {}
        for line in f:
            line_split = line.strip().split('\t')
            config_elem[line_split[0]] = line_split[1]
    # Step 1: Obtain names from ped file
    names = helper_functions.ped_file_reader(config_elem["PED_FILE"])

    # Step 2: Read in the vcf file
    df = helper_functions.read_VCF(config_elem["VCF_FILE"], names)

    # Step 3: obtain only SNPs and then group by chromosome name
    SNP_df = helper_functions.SNP_filter(df)
    snp_chr_dict = dict(tuple(SNP_df.groupby(["#CHROM"])))

    # Step 4: Generate results
    #start = time.time()
    mosaicism_outcome = [phasing_functions.runner(config_elem["PROBAND_NAME"], names[0], names[1],
                                                  int(config_elem["SAMPLE_SIZE"]), float(config_elem["T_THRES"]),
                                                  snp_chr_dict[chr_name], float(config_elem["EDGE_DETECTION_WIDTH"])) for chr_name in
                         range(1, 23)]
    #end = time.time()
    # print(" ".join(["It took", str((end - start) / 60), "minutes"]))
    # Step 5: Write results to file
    with open(config_elem["OUTPUT_FILE"], 'w') as output_file:
        output_file.write(" ".join(["Mosaicism results for", config_elem["PROBAND_NAME"]]))
        output_file.write("\n")
        output_file.write("\t".join(["chr_number", "VCF_start", "VCF_end"]))
        output_file.write("\n")
        for i, elem in enumerate(mosaicism_outcome):
            if elem is not None:
                output_file.write("\t".join([str(i + 1), str(elem[0]), str(elem[1])]))
                output_file.write("\n")


main('phasing_config_file.txt')

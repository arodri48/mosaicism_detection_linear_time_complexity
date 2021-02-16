from child import Child
import helper_functions


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

    # Step 3: obtain only SNPs
    SNP_df = helper_functions.SNP_filter(df)

    # New step 4: Analyze only the proband and feed in finder by chromosome
    chromosome_list = [i for i in range(1,23)]
    chromosome_list.append("MT")
    # TODO: make new constructor
    proband = Child(names[2], names[0], names[1], SNP_df)

    for chr_name in chromosome_list:
        # TODO: Write this function
        proband.runner(chr_name, config_elem["SAMPLE_SIZE"], config_elem["T_THRES"])

    # New Step 5: Once each chromosome is ran, check output
    if len(proband.results) == 0:
        print("No mosaicism detected")
    else:
        print("mosaicism detected")
        # TODO: Write output to file
        for elem in proband.results:
            print(elem)










    # Step 4: Create list of children
    children = [Child(names[i], names[0], names[1]) for i in range(2, len(names))]


    # Step 5: Generate phasable SNP data for each child
    for elem in children:
        elem.phasable_snp_determiner(chr7_df)

    # Step 6: Do sliding t-test for each child and see if possible mosaicism present
    for elem in children:
        elem.t_test_snps(10)

    # Step 7: Check if a child is mosaic
    mosaic_child_present = False
    for elem in children:
        if elem.index_diff_arr_start_of_mosaicism != 0:
            mosaic_child_present = True
            # detect region of mosaicism in the child

    # if all clear, print out
    if not mosaic_child_present:
        print("No child is mosaic")


main('phasing_config_file.txt')

from child import Child
import helper_functions


def main(config_file_path):
    # Step 0: Read config file
    with open(config_file_path, 'r') as f:
        config_elem = [line.strip().split('\t')[1] for line in f]
    # Step 1: Obtain names from ped file
    names = helper_functions.ped_file_reader(config_elem[0])

    # Step 2: Read in the vcf file
    df = helper_functions.read_VCF(config_elem[1], names)

    # Step 3: Filter by specific chromosome and obtain only SNPs
    chr7_df = helper_functions.filter_VCF_by_chr_and_SNP(df, 7)

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
        if elem.est_start_of_mosaicism != 0:
            mosaic_child_present = True
            # detect region of mosaicism in the child

    # if all clear, print out
    if not mosaic_child_present:
        print("No child is mosaic")


main('phasing_config_file.txt')

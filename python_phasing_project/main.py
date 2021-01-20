from child import Child
import helper_functions


def main():
    # Step 1: Obtain names from ped file
    names = helper_functions.ped_file_reader('UDP18111.ped')

    # Step 2: Read in the vcf file
    df = helper_functions.read_VCF('cohort_joint_genotyped_UDP18111.vcf', names)

    # Step 3: Filter by specific chromosome and obtain only SNPs
    chr7_df = helper_functions.filter_VCF_by_chr_and_SNP(df, 7)

    # Step 4: Create list of children
    children = [Child(names[i], names[0], names[1]) for i in range(2, len(names))]
    for elem in children:
        print(elem.name)

    # Step 5: Generate phasable SNP data for each child
    for elem in children:
        elem.phasable_snp_determiner(chr7_df)

    # Step 6: Do sliding t-test for each child and see if possible mosaicism present
    for elem in children:
        elem.t_test_snps(10)

    # Step 7: Check if a child is mosaic
    mosaic_child_present = False
    for elem in children:
        if (elem.est_start_of_mosaicism != 0):
            mosaic_child_present = True
            # detect region of mosaicism in the child

    # if all clear, print out
    if (not mosaic_child_present):
        print("No child is mosaic")


main()

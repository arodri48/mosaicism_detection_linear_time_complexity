package com.company;
import java.io.*;
import java.util.*;
import htsjdk.samtools.BAMFileReader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.SamLocusIterator;

public class PhasingProgram{
    public String [] patientNames;
    public ArrayList<String> proband_genotypes;
    public ArrayList<String> paternal_genotypes;
    public ArrayList<String> maternal_genotypes;

    //initialization of the counters
    public int headerCounterChr ; //initializes the variable vsLineReader to locate the chromosome column
    public int headerCounterMutType ;   //initializes the variable headerCounterMutType to count the number of headers to locate the column containing "muttype" data
    public int headerCounterVCF ;    //initializes the variable headerCounterVCF to count the number of headers to locate the column containing "VCF_Position" data
    public int headerCounterGenotype ;   //initializes the variable headerCounterGenotype to count the number of headers to locate the start of single patient data
    public int headerCounterTotalPeople ;   //initializes the variable headerCounterTotalPeople to count the total number of people in the file
    public int headerCounterColumns ; //initializes the variable counting the number of columns
    public int headerCounterColumns2 ; //initializes the second variable counting the number of columns
    public int patientCounter ; //initializes the variable that counts up the number of patients that will be added to the string array of all of the patients
    public int proband_index = 0;
    public int father_index = 0;
    public int mother_index = 0;

    public String [] header_info_finder(String curLine, String proband, String father, String mother) {//, int patientCounter, int headerCounterChr, int headerCounterMutType, int headerCounterVCF, int headerCounterGenotype, int headerCounterTotalPeople, int headerCounterColumns, int headerCounterColumns2
        String[] headers = curLine.split("\t"); //splits line by tab

        int columns = headers.length; //declaring an integer called "columns".  columns has a value that is equal to the total number of headers in the vs file.
        //this loop increases the header counter unless a header is named "Chr" is found or the end of the file has been reached to locate the first header containing single patient genotype data
        while (!(headers[headerCounterChr].equals("Chr"))&&!(headerCounterChr == columns)) {
            ++this.headerCounterChr;  //identifies the position of the "Chr" header column in the variant call file
        }
        //this loop increases the header counter unless a header is named "muttype" is found or the end of the file has been reached to locate the first header containing single patient genotype data
        while (!(headers[headerCounterMutType].equals("muttype"))&&!(headerCounterMutType == columns)) {
            ++this.headerCounterMutType;   //identifies the position of the "muttype" header column in the variant call file
        }
        //this loop increases the header counter unless a header is named "VCF_Position" is found or the end of the file has been reached to locate the first header containing single patient genotype data
        while (!(headers[headerCounterVCF].equals("VCF_Position"))&&!(headerCounterVCF == columns)) {
            ++this.headerCounterVCF;  //identifies the position of the "VCF_Position" header column in the variant call file
        }
        //this loop increases the header counter unless a header ending in .NA is found or the end of the file has been reached to locate the first header containing single patient genotype data
        while (!(headers[headerCounterGenotype].endsWith(".NA"))&&!(headerCounterGenotype == columns)) {
            ++this.headerCounterGenotype;     	//identifies the position of the first occurrence of single patient genotype data in the headers array
        }

        ArrayList<String> patientNames = new ArrayList<>();
        while(this.headerCounterColumns != columns) {
            if (headers[headerCounterColumns].endsWith(".NA")) {
                String patient_name = headers[headerCounterColumns].substring(0,headers[headerCounterColumns].length() - 3);
                patientNames.add(patient_name);
                if (patient_name.equals(proband)){
                    this.proband_index = headerCounterColumns;
                }
                else if(patient_name.equals(father)){
                    this.father_index = headerCounterColumns;
                }
                else if(patient_name.equals(mother)){
                    this.mother_index = headerCounterColumns;
                }
            }
            ++this.headerCounterColumns;
        }
        this.headerCounterTotalPeople = patientNames.size();
        this.patientCounter = patientNames.size();
        this.headerCounterColumns2 = this.headerCounterColumns;

        return patientNames.toArray(new String[0]);
    }

    public boolean het_check(String genotype_string){
        // first check if string is "NA"
        char first_letter = genotype_string.charAt(0);
        return  first_letter != genotype_string.charAt(1) && first_letter != 'N';
    }

    public boolean denovo_het_check(String proband_geno, String father_geno, String mother_geno){
        char first_letter = proband_geno.charAt(0);
        char second_letter = proband_geno.charAt(1);
        return (!father_geno.contains(""+first_letter) || !mother_geno.contains(""+second_letter)) && (!father_geno.contains(""+second_letter) || !mother_geno.contains(""+first_letter));
    }

    public Set<String> phasable_snp_finder(BufferedReader vsfile) throws IOException{
        String proband_genotype;
        String father_genotype;
        String mother_genotype;
        String curVCF_pos;
        String curLine;
        String [] split_line;
        String prevVCF = "00000";
        Set<String> duplicate_set = new HashSet<>();
        Set<String> phasable_snp_VCF_positions = new LinkedHashSet<>();
        while((curLine = vsfile.readLine()) != null){
            split_line = curLine.split("\t");
            proband_genotype = split_line[proband_index];
            father_genotype = split_line[father_index];
            mother_genotype = split_line[mother_index];
            curVCF_pos = split_line[headerCounterVCF];
            if (split_line[headerCounterMutType].equals("SNP")){
                // if variant is SNP, see if proband is a het
                if (het_check(proband_genotype)){
                    //if proband is het, first check if at least one of the parents homozygous (not phasable if both are het), but it is not the case that they are homozygous for the same allele
                    if (!(het_check(father_genotype) && het_check(mother_genotype)) && !father_genotype.equals(mother_genotype)){
                        // if none of the three are NA, save the VCF VCF_Position
                        if (!(proband_genotype.charAt(0) == 'N' || father_genotype.charAt(0) == 'N' || mother_genotype.charAt(0) == 'N')){
                            // check if proband alleles are in parents (that it's mendelian consistent)
                            if (!denovo_het_check(proband_genotype, father_genotype, mother_genotype)){
                                // check if prevVCF position matches the current one
                                if (!curVCF_pos.equals(prevVCF)){
                                    // first time seeing the position; add position to linkedhash set and then set prevVCF to new VCF position
                                    phasable_snp_VCF_positions.add(curVCF_pos);
                                    prevVCF = curVCF_pos;
                                }
                                else {
                                    // We have a duplicate SNP at the same position; add to duplicate set for removal from linked hash set later
                                    duplicate_set.add(curVCF_pos);
                                }
                            }


                            // 1 in a thousand, snp distance between each other on average
                        }
                    }
                }
            }
        }
        for (String value: duplicate_set){
            phasable_snp_VCF_positions.remove(value);
        }
        return phasable_snp_VCF_positions;
    }

    public void phasable_genotype_getter(BufferedReader input_file, Set<String> phasable_snps) throws IOException {
        ArrayList<String> proband_geno = new ArrayList<>(phasable_snps.size());
        ArrayList<String> paternal_geno = new ArrayList<>(phasable_snps.size());
        ArrayList<String> maternal_geno = new ArrayList<>(phasable_snps.size());
        String curLine;
        String [] split_line;
        String cur_phasable_pos;
        Iterator itr = phasable_snps.iterator();
        curLine = input_file.readLine();
        split_line = curLine.split("\t");
        cur_phasable_pos = (String) itr.next();
        long int_cur_pos;
        long pos_vs_file;
        while (true){
            int_cur_pos = Long.valueOf(cur_phasable_pos);
            pos_vs_file = Long.valueOf(split_line[headerCounterVCF]);
            if (int_cur_pos > pos_vs_file){
                curLine = input_file.readLine();
                if (curLine != null){
                    split_line = curLine.split("\t");
                }
                else {
                    break;
                }
            }
            else if (int_cur_pos == pos_vs_file){
                // add to the genotypes
                proband_geno.add(split_line[proband_index]);
                paternal_geno.add(split_line[father_index]);
                maternal_geno.add(split_line[mother_index]);
                // then check
                curLine = input_file.readLine();
                if (curLine != null && itr.hasNext()){
                    cur_phasable_pos = (String) itr.next();
                    split_line = curLine.split("\t");
                }
                else {
                    break;
                }
            }
            else {
                if (itr.hasNext()){
                    cur_phasable_pos = (String) itr.next();
                }
                else {
                    break;
                }
            }
        }
        this.proband_genotypes = proband_geno;
        this.paternal_genotypes = paternal_geno;
        this.maternal_genotypes = maternal_geno;

    }

    public char [][] phasing_determiner(ArrayList<String> proband, ArrayList<String> father, ArrayList<String> mother){
        // Step 1: Initialize N x 2 char array; first column is father allele; second is mother allele
        char [][] phased_array = new char[proband.size()][2];
        // Step 2: Iterate through the arraylists and fill in the array
        int snp_len = proband.size();
        char first_allele;
        char second_allele;
        String proband_geno;
        String father_geno;
        String mother_geno;
        for (int i = 0; i != snp_len; ++i){
            proband_geno = proband.get(i);
            father_geno = father.get(i);
            mother_geno = mother.get(i);
            first_allele = proband_geno.charAt(0);
            second_allele = proband_geno.charAt(1);
            int father_first_allele_indexof = father_geno.indexOf(first_allele);
            // first check if first allele letter is in both parents; if so, we are dealing with f0xf1 cross
            if (father_first_allele_indexof != -1 && mother_geno.indexOf(first_allele) != -1){
                // check if second allele is in father
                if (father_geno.indexOf(second_allele) != -1){
                    // second allele is in father, so first allele belongs to mother
                    phased_array[i][0] = second_allele;
                    phased_array[i][1] = first_allele;
                }
                else {
                    phased_array[i][0] = first_allele;
                    phased_array[i][1] = second_allele;
                }
            }
            // first allele is in only one of the parents
            else {
                if (father_first_allele_indexof != -1){
                    phased_array[i][0] = first_allele;
                    phased_array[i][1] = second_allele;
                }
                else {
                    phased_array[i][0] = second_allele;
                    phased_array[i][1] = first_allele;
                }
            }
        }
        return phased_array;
    }

    public int [][] read_depth_finder(Set<String> snp_vcf_pos, char [][] allele_parents, String father_bam, String mother_bam){
        int len_SNPS = snp_vcf_pos.size();
        int [][] read_depths_array = new int[len_SNPS][2];
        // First get the read depths for the father
        // Step 1: Create iterators for SNP_VCF_POS and the father's BAM
        Iterator it1 = snp_vcf_pos.iterator();
        SamReader reader = SamReaderFactory.makeDefault().open(new File(father_bam));
        SamLocusIterator sli1 = new SamLocusIterator(reader);
        // Step 2: Load the first values of the iterators
        String cur_phasable_pos = (String) it1.next();
        SamLocusIterator.LocusInfo li = sli1.next();
        long int_cur_pos;
        long pos_bam_file;
        int counter1 = 0;
        while (true){
            int_cur_pos = Long.valueOf(cur_phasable_pos);
            pos_bam_file = li.getPosition();
            if (int_cur_pos > pos_bam_file){
                if (sli1.hasNext()){
                    li = sli1.next();
                }
                else {
                    break;
                }
            }
            else if (int_cur_pos == pos_bam_file){
                // TODO: if positions match, check if line is proper allele; if so, save the read depth and increment
                //  counter; if not proper allele, read the next line in the bam file
            }
            else {
                if (it1.hasNext()){
                    cur_phasable_pos = (String) it1.next();
                }
                else {
                    break;
                }
            }
        }
        reader.close();


        // Second, get the read depths for the mother
        return read_depths_array;
    }

    public void runner(String VS_file, String BAM_file, int chunk_size, String proband, String father, String mother, String father_bam, String mother_bam) throws IOException {
        // Step 1: Read in VS file and determine where proband is a het at
        BufferedReader vsfile_reader = new BufferedReader(new FileReader(VS_file));
        String curLine;
        curLine = vsfile_reader.readLine();
        this.patientNames = header_info_finder(curLine,proband,father, mother);
        Set<String> phasable_snp_positions = phasable_snp_finder(vsfile_reader);
        vsfile_reader.close();

        // Step 2: Read in the VS file again, this time to obtain the genotypes at the right positions
        BufferedReader vs2_reader =  new BufferedReader(new FileReader(VS_file));
        curLine = vs2_reader.readLine();
        phasable_genotype_getter(vs2_reader, phasable_snp_positions);
        vs2_reader.close();

        // Step 3: Now that we have the positions and genotypes at each location, it's time to do the phasing; create
        // N x 2 String Array, where the first column is the allele from the dad and the second is the allele from
        char [][] parent_allele_array = phasing_determiner(proband_genotypes, paternal_genotypes, maternal_genotypes);

        // Step 4: Obtain the read depth of the allele from the father and the mother from the BAM files and put into
        // N x 2 array
        int [][] read_depth_array = read_depth_finder(phasable_snp_positions, parent_allele_array, father_bam, mother_bam);







        // Step 2: Samtools to get read depths
        // open up the BAM file
        SamReader reader = SamReaderFactory.makeDefault().open(new File(BAM_file));
        //create an iterator
        SamLocusIterator sli = new SamLocusIterator(reader);
        int bam_position;
        // for each phasable, valid SNP position
        for (String val : phasable_snp_positions) {
            // for each line in the BAM file
            for (SamLocusIterator.LocusInfo li : sli) {
                // check if position is same
                bam_position = li.getPosition();
                if (bam_position == val){
                    // it's the same position
                }
                else if (bam_position < val) {
                    // bam position is less than current phasable position
                }
                else {
                    // bam position is greater than current phasable position
                }

            }
        }
        // Step 3: Gang maternal and paternal read depths;

    }
    public static void main(String []args ) throws IOException{
        PhasingProgram obj = new PhasingProgram();
        BufferedReader config_file = new BufferedReader(new FileReader(args[0]));
        String curLine;
        String [] split_line;
        String [] values = new String[6];
        for (int i = 0; i != 6; ++i){
            curLine = config_file.readLine();
            split_line = curLine.split("\t");
            values[i] = split_line[1];
        }
        config_file.close();
        //obj.runner(values[0], values[1], Integer.parseInt(values[2]), values[3], values[4], values[5]);
        System.gc();
        System.exit(0);
    }
}

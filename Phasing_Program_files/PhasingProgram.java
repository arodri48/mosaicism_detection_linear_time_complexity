import java.io.*;
import java.util.*;
class PhasingProgram{
  String [] patientNames;
	Map <String, Integer> patientname_indices;

	//initialization of the counters
	public int headerCounterChr ; //initializes the variable vsLineReader to locate the chromosome column
	public int headerCounterMutType ;   //initializes the variable headerCounterMutType to count the number of headers to locate the column containing "muttype" data
	public int headerCounterVCF ;    //initializes the variable headerCounterVCF to count the number of headers to locate the column containing "VCF_Position" data
	public int headerCounterGenotype ;   //initializes the variable headerCounterGenotype to count the number of headers to locate the start of single patient data
	public int headerCounterTotalPeople ;   //initializes the variable headerCounterTotalPeople to count the total number of people in the file
	public int headerCounterColumns ; //initializes the variable counting the number of columns
	public int headerCounterColumns2 ; //initializes the second variable counting the number of columns
	public int patientCounter ; //initializes the variable that counts up the number of patients that will be added to the string array of all of the patients
	public int configFileCounter = 0; //counts up the number of lines in the  file
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

		ArrayList<String> patientNames = new ArrayList<String>();
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
		if (!(first_letter == 'N' || first_letter == genotype_string.charAt(1))){
			return true;
		}
		return false;
	}

  public void runner(String VS_file, String BAM_file, int chunk_size, String proband, String father, String mother) throws IOException {
    // Step 1: Read in VS file and determine where proband is a het at
    BufferedReader vsfile_reader = new BufferedReader(new FileReader(VS_file));
    String curLine;
    String [] split_line;
    curLine = vsfile_reader.readLine();
    this.patientNames = header_info_finder(curLine, proband, father, mother);
    while((curLine = vsfile_reader.readLine()) != null){
      split_line = curLine.split("\t");
      // check if variant is snp
			if (split_line[headerCounterMutType].equals("SNP")){
				// if variant is SNP, see if proband is a het
				if (het_check(split_line[proband_index])){
					//if proband is het, first check if at least one of the parents homozygous (not phasable if both are het)
					if (!(het_check(split_line[father_index]) && het_check(split_line[mother_index]))){
						// if none of the three are NA, save the VCF VCF_Position
						if (!(split_line[proband_index].charAt(0) == 'N' || split_line[father_index].charAt(0) == 'N' || split_line[mother_index].charAt(0) == 'N')){
							// save the VCF position

						}
					}
				}
			}

    }
		vsfile_reader.close();
  }


  public static void main(String args[]) throws IOException{
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
    obj.runner(values[0], values[1], Integer.parseInt(values[2]), values[3], values[4], values[5]);
    System.gc();
    System.exit(0);
  }
}

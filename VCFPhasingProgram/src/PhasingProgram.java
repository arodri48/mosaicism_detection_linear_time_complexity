import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

public class PhasingProgram {
	public int proband_index = 0;
	public int father_index = 0;
	public int mother_index = 0;
	public int chrom_index = 0;
	public int pos_index = 0;
	public int qual_index = 0;
	public int ref_index = 0;
	public int alt_index = 0;
	public ArrayList<Integer> paternal_read_depth;
	public ArrayList<Integer> maternal_read_depth;

	public void header_info_finder(String header_line, String proband, String father, String mother){
		String [] headers = header_line.split("\t");
		int col_num = headers.length;
		// find index of '#CHROM'
		while(!(headers[this.chrom_index].equals("#CHROM")) && this.chrom_index != col_num){
			++this.chrom_index;
		}
		// find index of 'POS'
		while(!(headers[this.pos_index].equals("POS")) && this.pos_index != col_num){
			++this.pos_index;
		}
		// find index of 'REF'
		while(!(headers[this.ref_index].equals("REF")) && this.ref_index != col_num){
			++this.ref_index;
		}
		// find index of 'ALT'
		while(!(headers[this.alt_index].equals("ALT")) && this.alt_index != col_num){
			++this.alt_index;
		}
		// find index of 'QUAL'
		while(!(headers[this.qual_index].equals("QUAL")) && this.qual_index != col_num){
			++this.qual_index;
		}
		// find index of proband
		while(!(headers[this.proband_index].equals(proband)) && this.proband_index != col_num){
			++proband_index;
		}
		// find index of father
		while(!(headers[this.father_index].equals(father)) && this.father_index != col_num){
			++father_index;
		}
		// find index of mother
		while(!(headers[this.mother_index].equals(mother)) && this.mother_index != col_num){
			++mother_index;
		}

	}
	public String [] ped_name_extractor(String ped_file) throws IOException{
		String [] names = new String[3];
		String curLine;
		String [] split_line;
		BufferedReader ped_file_reader = new BufferedReader(new FileReader(ped_file));
		curLine = ped_file_reader.readLine();
		split_line = curLine.split("\t");
		// first entry in names is proband, second is father, third is mother
		names[0] = split_line[0];
		names[1] = split_line[2];
		names[2] = split_line[3];
		ped_file_reader.close();
		return names;
	}
	public void read_depth_determiner(BufferedReader VCF_file_reader){

	}
	public void runner(String VCF_file_path, String ped_file, String output_path) throws IOException {
		// Step 1: Open ped file and extract proband, father, and mother names
		String [] names = ped_name_extractor(ped_file);
		// Step 2: Open the VCF file and extract the names and other information from the header
		BufferedReader VCF_file_reader = new BufferedReader(new FileReader(VCF_file_path));
		String curLine;
		while ((curLine = VCF_file_reader.readLine()) != null){
			if (curLine.charAt(1) != '#'){
				break;
			}
		}
		header_info_finder(curLine, names[0], names[1], names[2]);
		// Step 3: At each location that satisfies criteria for child and parents, extract information about read depths
		// for alleles
		read_depth_determiner(VCF_file_reader);

		// Step 4: Calculate super SNPs and save to file
	}
	public static void main(String [] args) throws IOException{
		// Step 1: Read parameters from config file
		BufferedReader config_file = new BufferedReader(new FileReader(args[0]));
		String curLine;
		String [] split_line;
		String [] values = new String[5];
		for (int i = 0; i != 5; ++i){
			curLine = config_file.readLine();
			split_line = curLine.split("\t");
			values[i] = split_line[1];
		}
		config_file.close();


		PhasingProgram snp_finder = new PhasingProgram();
		snp_finder.runner("Path", "Path", "Path");
		System.gc();
		System.exit(0);
	}
}

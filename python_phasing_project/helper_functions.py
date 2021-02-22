import io
import pandas as pd


def read_VCF(vcf_path, name_list):
    with open(vcf_path, 'r') as f:
        lines = [line for line in f if not line.startswith('##')]
    df = pd.read_csv(io.StringIO(''.join(lines)), sep='\t')
    df = df[~df['REF'].str.contains(",")]
    df = df[~df['ALT'].str.contains(",")]
    df['REF'] = df['REF'].astype('str')
    df['ALT'] = df['ALT'].astype('str')
    for elem in name_list:
        df[elem] = df[elem].astype('str')
    df['QUAL'] = pd.to_numeric(df['QUAL'], errors='coerce')
    df['POS'] = pd.to_numeric(df['POS'], errors='coerce')
    return df


def filter_VCF_by_chr_and_SNP(df, chr_number):
    return df[(df["#CHROM"] == chr_number) & (df['REF'].str.len() == 1) & (df['ALT'].str.len() == 1) & (
            df['QUAL'] > 100)].reset_index(drop=True)


def SNP_filter(df):
    return df[(df['REF'].str.len() == 1) & (df['ALT'].str.len() == 1) & (
            df['QUAL'] > 100)].reset_index(drop=True)


def chromosome_filter(df, chr_name):
    return df[(df["#CHROM"] == chr_name)].reset_index(drop=True)


def ped_file_reader(ped_file_path):
    with open(ped_file_path, 'r') as f:
        names = [l.strip().split('\t')[1] for l in f]
    actual_names = [names[1], names[2], names[0]]
    for i in range(3, len(names)):
        actual_names.append(names[i])
    return actual_names

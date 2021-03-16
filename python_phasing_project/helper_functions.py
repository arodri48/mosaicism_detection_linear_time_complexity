import io
import pandas as pd


def read_VCF(vcf_path, name_list):
    with open(vcf_path, 'r') as f:
        lines = [line for line in f if not line.startswith('##')]
    dtype_dic = {'REF': 'str', 'ALT': 'str'}
    for elem in name_list:
        dtype_dic[elem] = 'str'
    df = pd.read_csv(io.StringIO(''.join(lines)), sep='\t', dtype=dtype_dic)
    df['QUAL'] = pd.to_numeric(df['QUAL'], errors='coerce')
    df['POS'] = pd.to_numeric(df['POS'], errors='coerce')
    return df


def SNP_filter(df):
    return df[(df['REF'].str.len() == 1) & (df['ALT'].str.len() == 1) & (
            df['QUAL'] > 100)]


def ped_file_reader(ped_file_path):
    with open(ped_file_path, 'r') as f:
        names = [l.strip().split('\t')[1] for l in f]
    actual_names = [names[1], names[2], names[0]]
    for i in range(3, len(names)):
        actual_names.append(names[i])
    return actual_names

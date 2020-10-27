import argparse
import matplotlib.pyplot as plt
import pandas as pd
import os
def main():
	# read in arguments
	parser = argparse.ArgumentParser(prefix_chars='-')
	parser.add_argument('-i', required=True)
	parser.add_argument('-o', required=True)
	args = parser.parse_args()
	# load in the data
	data = pd.read_csv(args.i, sep='\t', header=None)
	# create plot
	rows = data[0].count()
	row_indexes = [i for i in range(rows)]
	plt.scatter(row_indexes, data[0], s=1)
	plt.ylim(0.15, 0.85)
	plt.xlabel('SNP')
	plt.ylabel('BAF')
	# save to file
	basename = os.path.splitext(os.path.basename(args.i))[0]
	image_path = args.o + basename +".png"
	plt.savefig(image_path, dpi=300)

main()

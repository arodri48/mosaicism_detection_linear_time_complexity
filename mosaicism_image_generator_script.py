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
	data = pd.read_csv(args.i, sep='\t')
	# create plot
	plt.plot(data['Start_Pos'], data['Discriminant_Value'])
	plt.ylim(0.0, 9)
	plt.xlabel('Start Position')
	plt.ylabel('Discriminant Value')
	# save to file
	basename = os.path.splitext(os.path.basename(args.i))[0]
	image_path = args.o + basename +".png"
	plt.savefig(image_path, dpi=300)

main()

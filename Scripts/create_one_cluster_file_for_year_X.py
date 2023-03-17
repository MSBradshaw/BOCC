import os
import argparse

# create a function to parse the command line arguments
def parse_args() -> argparse.Namespace:
	parser = argparse.ArgumentParser(description='A simple command line interface for the Bocc library')
	parser.add_argument('-y', '--year', required=True, help='Year to process')
	return parser.parse_args()

args=parse_args()

with open(f'all_paris_coms.{args.year}.txt','w') as outfile:
	for file in os.listdir(f'SubComs/{args.year}/'):
		if 'paris' not in file:
			continue
		print(file)
		prefix = file.split('.coms.')[0]
		for line in open(f'SubComs/{args.year}/'+file,'r'):
			outfile.write(prefix + ':' + line)

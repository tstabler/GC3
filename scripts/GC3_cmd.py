#!/usr/local/packages/python-3.8.2/bin/python3.8

from GC3 import *
import argparse
import os

###########################################################

print('START')

parser = argparse.ArgumentParser(description='Process .gz file.')
parser.add_argument('-v', '--ascii-file',type=str,required=True,help='Path to input file')
parser.add_argument('-o', '--output', help = 'Directs the output to directory of your choice')
parser.add_argument('-s', '--start',required=True, help = 'User-defined start coorindate')
parser.add_argument('-e', '--end',required=True, help = 'User-defined end coordinate')
parser.add_argument('-m', '--molecule',required=True, help = 'User-defined molecule')
parser.add_argument('-i', '--interval', help = 'Interval size')
parser.add_argument('-j', '--step', help = 'Step size')
parser.add_argument('-f', '--function',required=True, help = 'GC3 function option')

arg = parser.parse_args()

file_name=os.path.basename(arg.ascii_file)


#Use GC3 functions to process coverage data - User parameters required between < and >. Do not include < >!
with open(arg.ascii_file) as file:
	print("Starting function...")
	if arg.function=="window":
		GC3_output = calc_coverage_int(file, start_position=int(arg.start), end_position=int(arg.end), molecule=arg.molecule, interval=int(arg.interval), jump=int(arg.step))
	elif arg.function=="mean":
		GC3_output = mean_coverage(file, start_position=int(arg.start), end_position=int(arg.end), molecule=arg.molecule)

with open(arg.output, 'a') as output1:
	output1.write(f"{file_name}, {GC3_output} \n")
	output1.close()

print("END")

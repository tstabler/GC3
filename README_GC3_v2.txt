README FILE FOR GC3 v2.0

SEP2022
Thomas C. Stabler

#############################################################################################

ABOUT

This file contains instructions regarding the implementation and operation of the process termed "Gene Coverage and Classification" or GC3, a locus coverage assessment tool. There are two junctions in the GC3 process that requires user defined parameters. Steps are detailed below with example commands. 

All text between < and > refers to where user input is required. < and > should NOT be included when scripts are run. 

#############################################################################################

PROGRAM DEPENDENCIES

Python v3 and above

R v4.1 and above
	Libraries: tidyr, ggplot2, writexl, dplyr, tibble, readxl, and reshape2

#############################################################################################

REQUIRED FILES

BED file containing coverage per position

To see how BED coverage files were generated, visit: 
https://github.com/igs-jcsilva-lab/variant-calling-pipelines

or also acceptable

a tab-delimited with with columns in the following order: molecule identifier (e.g. Pf3D7_08_v3),
	chromosomal position, and coverage value

#############################################################################################

To start, create a sample directory, followed by country, unique sample ID, and alignments.
Store all coverage files in their appropriate location according to this organization. 

NOTE: "samples" and "alignment" in the directory pathway are used as keywords to orientate GC3. To change keywords requires
amendment of GC3_delshell.sh script 

/samples/<country>/<sample_ID>/alignments/<coverage_file_name> 

Example:

#Making directory
mkdir -p /p_falciparum/samples/cambodia/IGS-CBD-001/alignments/

#############################################################################################

INSTRUCTIONS FOR GC3 PYTHON SCRIPT

The driver of GC3 is the shell command script (GC3_delshell.sh) that imports GC3 command script (GC3_cmd.py) and functions (GC3.py), 
parses coverage data and applies functions as defined by user. For analysis of coverage data that doesn't require significant memory 
or computational power, the command line in the shell command script can be directly applied into the command line. 

To adjust any parameters, user only needs to edit the command shell script and rerun it.

The following parameters can be input into GC3 (required inputs are marked with an "R"):

	-v (R) input_file (either direct path to coverage file or .txt file with list of pathways to coverage file(s))
	-o output_file
	-s (R) start coordinate
	-e (R) end coordinate
	-m (R) molecule 
	-i interval size
	-j step size
	-f (R) function -> "window" or "mean"


Examples are below for a commands used for sliding window output, positional coverage output, and mean coverage
output. These commands are either present in the command shell script or on the command line. 

It is recommended the GC3_delshell.sh script be run on the grid using a qsub command, but can be run locally. 
Make sure pathway to python package is correct in the GC3_delshell.sh script


If running directly on command line, input file needs to be unzipped (does not end in .gz) and pathway to python package defined. 

Example(s):

gzip -d <pathway/to/file/position-coverage.gz>

py=<"/usr/packages/python-3.5/bin/python-3.5">

Sliding window:

	$py GC3_cmd.py -v position-coverage -s 2838727 -e 2843703 -m "Pf3D7_13_v3" -i 500 -j 250 -f "window" -o output.txt
	
Positional coverage (extract coverage at every position - j input not required):

	$py GC3_cmd.py -v position-coverage -s 2838727 -e 2843703 -m "Pf3D7_13_v3" -i 1 -f "window" -o output.txt

	
Mean coverage between start and stop coordinates:

	$py GC3_cmd.py -v position-coverage -s 2838727 -e 2843703 -m "Pf3D7_13_v3" -f "mean" -o output.txt
	


Above commands are suitable if running a single coverage file through GC3. If running multiple coverage 
files through GC3, the best way is to generate a .txt file where each line is a directory to one a 
coverage file to be analyzed (Example .txt file included on https://github.com/igs-jcsilva-lab. 
User should use the GC3_delshell.sh script which includes a "for loop" through each line of the .txt file.

If using the shell script, the following variables will need to be defined:

	folder="<folder_where_txt_file_of_samples_is_located>"
	path="<pathway_to_folder>"
	file="<name_of_txt_file>" #Don't include .txt at the end of file name
	tag="<name_of_folder_where_output_files_go>"
 
Although this adds another layer, this ultimately makes it easier to use when analyzing coverage of multiple samples since user
will not need to define pathways to packages each time GC3 is run.  

Example of "For Loop" in shell script:

py=<"/usr/packages/python-3.5/bin/python-3.5">

samples=$(while read p; do echo "$p"; done < <txt_file_with_pathways_to_file>)

for i in $samples; 
do 
	echo "Starting file processing..."
	label=$(echo $i | sed 's/\//_/g' | sed -e 's/.*samples_\(.*\)_alignments.*/\1/')
	echo $label
	zcat $i | awk '$1=="Pf3D7_08_v3" && $2>=1364967 && $2<=1384877' > "$path"/"$tag"/"$label"_tmp
	$py GC3_cmd.py -v "$path"/"$tag"/"$label"_tmp -s 1372236 -e 1377299 -m "Pf3D7_08_v3" -i 1 -j 1 -f "window" -o "$path"/"$tag"/output.txt
	rm "$path"/"$tag"/"$label"_tmp
	echo "Completed samples $label"
done


There are several built in commands built into GC3 and GC3_cmd.py scripts that report progress of GC3. 
User should see "END" or "END SHELL SCRIPT" message at the end of any output error file. There should be only one output file per
run of GC3. This intermediate file will then become the input into the R script (Next junction)
#############################################################################################

INSTRUCTIONS FOR R SCRIPT

R script instructions are specific to being run locally on a computer, but can be run on the grid with amendments.
 
NOTE: If using R for the first time, Users will need to uncomment (remove "#" at the start of the line)
under the "Introduction" section below "Import libraries" all lines with "install.packages" command. 
This will install the required R libraries. Once it is run the first time, these commands can be commented out (insert "#" at start of each line) as before. 

Users will only need to provide input under the "User input" section of the R script. Instructions are
included within the R script as comments, but are also below. 

directory="<C:/pathway/to/intermediate/file>"
analysis="<C:/pathway/to/where/analysis/outputs/go>"

larger_input="<name_of_file_for_sliding_window_function>"

smaller_input="<name_of_file_if_function_interval_is_1>"

standardization="<name_of_file_for_mean_coverage_function>"

gene_start=<start_coordinate>
gene_end=<end_coordinate>

NOTE: if gene does not have an intron, set intron start/end to zero (0)
intron_start=<intron_start>
intron_end=<intron_end>

######
Optional parameters (can only be used if input includes targetted positions):

#Flanking genes
Set to 0 if not needed. 
NOTE: make sure coordinates are within intermediate file from GC3 python script


#Zoom
Set to 0 if not needed
NOTE: this is to look at areas of interest. For example, a gene's exon. 


#Subgrouping sample set
NOTE: R script currently set up to analyze 2 subgroups. Relevant lines will need to be uncommented to use.
If not needed, user will need to comment out lines again. 
Slightly tedious, but user will need to list samples that make up one of the subgroups

Format of list of samples need to be <country_sampleID>. User will then define subgroup names.

sample_change1=c("<country_sampleID_1>", "<country_sampleID_2>", etc.)

subgroup1 <- "Name_of_Group_1"
subgroup2 <- "Name_of_Group_2"

NOTE: when setting up subgroups, plots visuals can be weird. User can define y-limit values of plots 
to determine the best plot perspective. There is a y-limit values for non-standardized plots (limity)
and standardized plots (limity_stan)  


Once all necessary parameters are set, user will select entire R script (ctrl a) and run.
#############################################################################################

Outputs will be generated and GC3 is complete! 






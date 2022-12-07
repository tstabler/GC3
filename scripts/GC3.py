#!/usr/local/packages/python-3.8.5/bin/python3.8

#written by Tom Stabler, Msc (SwissTPH)


def calc_coverage_int(file, start_position=1, end_position=1000, molecule="Pf3D7_08_v3", interval=500, jump=100):
	"""Calculate average gene coverage from start_position to end_position by interval X every Y bp. Default is to do this on Chromosome 8"""
	coverage_list = {} #Create empty dictionary for positions and coverage
	tmp_start = start_position
	tmp_end = start_position + interval
	coverage_interval = []
	print("Start line read")

	if interval>1:
		while True:
			if (tmp_end > end_position):
				print("End position reached!")
				break
			
			for line in file:
				if not line.startswith("#"): #If headers inserted at first row
					linelist = line.split()
					Chrom, Pos, Coverage = linelist[0:3] #Assign header names to each line
					
					if Chrom == molecule: 
				
						if (tmp_start <= int(Pos) and int(Pos) < tmp_end): #If within interval window then parse out format to get DP and other info (still need to add GT)
							DP = int(float(Coverage)) 
							coverage_interval.append(DP)
							continue
								
						if (int(Pos) < tmp_start or int(Pos) >= tmp_end): 
							if (int(Pos) < tmp_start):
								continue 
								
							if (int(Pos) >= tmp_end):
								coverage = ((sum(coverage_interval)) / (interval))
								coverage_list[f'{tmp_start}:{tmp_end}'] = coverage
								coverage_interval.clear()
								tmp_start += jump
								tmp_end += jump
								file.seek(0) #Return to beginning of file 
								break
					
			else:
				print("Goodbye!")
				file.seek(0)
				break
	
	
	if interval==1:
		while True:
			for line in file:
				if not line.startswith("#"): #If headers assigned in file
					linelist = line.split()
					Chrom, Pos, Coverage = linelist[0:3] #Assign header names to each line
					
					if Chrom == molecule: #Restart for loop when scanning over entire genome
				
						if start_position <= int(Pos) <= end_position:
							coverage_list[f'{Pos}'] = int(Coverage)
							#print(f"For {Pos} coverage = {Coverage}")
							continue

						if int(Pos) > end_position: 
							#file.seek(0) #Return to beginning of file 
							print("End position reached!")
							break
					
			break
	print("Sliding window complete!")
	file.seek(0) #Return to beginning of file
	return coverage_list

def mean_coverage(file, start=1, end=1000000, molecule="Pf3D7_08_v3"):
	"""Calculate mean coverage between start and end positions"""
	coverage_list = {} #Create empty dictionary for positions and coverage
	coverage_interval = []
	for line in file:
		if not line.startswith("#"): #If headers inserted at first row
			linelist = line.split()
			Chrom, Pos, Coverage = linelist[0:3] #Assign header names to each line
			
			if Chrom == molecule:
		
				if (start <= int(Pos) and int(Pos) <= end): #If within start/end positions then parse out format to get DP and other info
					DP = int(float(Coverage)) 
					coverage_interval.append(DP)
					continue
						
				if (int(Pos) < start or int(Pos) >= end): 
					if (int(Pos) < start):
						continue 
						
					if (int(Pos) > end):
						#file.seek(0) #Return to beginning of file
						print("End of desired positions reached!")
						break
	coverage = ((sum(coverage_interval)) / ((end-start)+1))
	coverage_list[f'{start}:{end}'] = coverage
	print("Calculated mean coverage")
	print("Goodbye!")
	file.seek(0) #Return to beginning of file
	return coverage_list


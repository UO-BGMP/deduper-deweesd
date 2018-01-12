#!/usr/bin/env python3

################################################## Argparse setup ###########################################################

import argparse
import re


parser = argparse.ArgumentParser(description= 'Iterates through sorted SAM file and removes PCR duplicates from uniq_mapped reads (single & paired end). quality score per base was used to separate output file while UMIs vs known UMIs were also used to separate output file. Randomers were ignored)

parser.add_argument('-f', '--file_path', help= "specifying file path for SAM file and read in, r", required=True, type=str)
parser.add_argument('-u', '--UMI', help= "Optional flag designating file containing accepted UMIs from known UMIs list", type=str)
parser.add_argument('-p', '--paired_end', help=" File is paired end (not single_end)", required=False, action='store_true, type=str)
parser.add_argument('-h', '--help', help= "prints informative help message about script/output", required=True, type=str)
args = parser.parse_args()


# When known_UMIs is added #

if args.UMI:
	known_UMIs = []
	with open(args.UMI, 'r+') as umi_document:
		
		# for each line in known umi_document #
		for line in umi_document:
			
			# append each UMI to empty list "known_UMIs" #
			known_UMIs.appened(line.strip())

# When known_UMIs is not given #
else:
	list_not_given = True




# If paired_end option is exercised, exit run with error #
if args.paired_end:
    raise NameError('Support for paired_end reads is not valid @ this time')

################################# logical grouping of data and functions (methods) via 'class' ######################################

class sort_sam_file():
    '''Using class & methods to iterate and filter through SAM file. 
    Process includes:
        
        - monitors & adjusts alignment start postion (via soft clipping amount ~ CIGAR string "S").
        - converts quality into phred score (Using ASCII table ~ assuming ASCII being 33).
        - checks for strandeness and unique mapping.
        - checks chromosome positioning
        - grabs molecular id (umi barcode).
    
    
    SAM file have the following properties:
        
        Attributes:
            line: looking at current_line (single) from sam file (string)
            start_position: left most mapping postion on said line ~ soft clip adjusted (integer)
            convert_phred: line quality score (integer)
            strand (Bit-wise flag): lines strandeness (string)
            soft_clipped: determines if seq featured soft clipping (integer)
            chromosome: which chromosome (integer)
    '''
    
    
    ## using '_init_' to construct the class (initializing atrributes listed above) and 'self' to call the instance of the object within class ##
    def __init__(self, current_line):
        
        self.current_line = line()
        self.start_position = int(line.split('\t')[3])
        self.q_score = int(line.split('\t')[4])
        self.bit_flag = int(line.split('\t')[1])
        self.soft_clipped = line.split('\t')[5]
        self.chromosome = int(line.split('\t')[2])
        


## function to check for soft_clipping ##    
def soft_checker(start_position):
    """ function to check for soft clipping *returns Soft clip amt"""
    
    ## calling cigar character "S" and split in-line ##
    cigar_string = cigar_string.split("S")
    
    ## using re to filter through all possible characters and cigar string as first character, ignore ##
    if re.search('[A-Z]', cigar_string[0]):
        pass
        
    ## taking start position as the int value minus the value prior to cigar letter (i.e., "S") ~ method to remove soft clippings in file ##
    else:
        start_position = int(start_position) - int(cigar_string[0])
            
    return start_position



def umi_checker(UMI):
    '''cross-reference current_UMI to known_UMI list ~ site: https://github.com/Leslie-C/Deduper/blob/master/STL96.txt'''
        
    if UMI in known_UMIs:
        return True
    
    else:
        # proceeds to ignore and throw this read out of file # 
        return False



def bit_checker(bit_flag, paired_end):
    '''Iterates through bitwise flag to check for strandedness (+/-).  
    Assumes reads are uniquely mapped, otherwise returns "None". Assumes data are single-end. Returns 
    "+" or "-", depending on strand.'''
        
    ## if bit flag 0x4 is not present, read is mapped ##
    if ((bit_flag & 4) != 4):
        mapped = True
    else:
    	#raise NameError("read is unmapped")
        return None 

## 0x10 is SEQ of reverse complemented (negative strand if TRUE) ##
    strand = "+"
    if (bit_flag & 16) == 16: # if bitwise flag 0x10 is true, strand is negative (reverse)
        strand == "-"
    return strand
    
    
## Using bitwise flag to determine if reads are paired_ended. Returns the corresponding statements ("forward/reverse" strands).    
    
    if paired_end == True:
    	 
    	#Checking if bit 64 (first segment) is "True"
    	if (bit_flag & 40) == 40:
    		current_read = "forward_strand"
    	
    	# checking if bit 128 (last segment) is "True"
    	if (bit_flag & 80) == 80:
    		current_read = "reverse_strand"
		
		# If true, return current_read
    	return current_read


####################################################### Setup for iteration & outputfil format ############################################################

## Set counter to 0 in order to tally up all duplicates from file as well as count for un-used UMIs ## 
duplicates_cnt = 0



## Generate output_sam file w.out duplications in said file ## 
duplicated = open(pre_fx+"_deduped.sam",'w')

## calling the prefix in file to change name output ##
pre_fx = args.filename[:-4]


## set file to absolute path ##
file = "/Users/Deweesd/Desktop/DeDuper_PCR/test_1.sam"


####################################################### Iterate through SAM file ############################################################

with open(file, 'r') as sam_file_1:
    # setting line count to 0 #
    line_cnt = 0
    
    ## looking at each line in SAM file (fh) ##
    for read_line in sam_file_1:
        
        ## iterating through SAM file (HEADER) to see if line starts with certain character (@), if so, ignore and move onto next line##
        if read_line.startswith('@'):
            duplicated.write(read_line)
            continue
        
        elif read_line.startswith('@') == False:
        	# for each line, increment by 1 #
        	line_cnt += 1
        	
        	## looking at each line that post 'header label' and strips by new line and splits by tab ##
        	current_line = read_line.strip('\n').split("\t")

        	
        	
        	############# Initiate reference variables for iteration ##############
        	
        	## Unadjusted raw_pos - soft clipped not checked ##
            pos = int(line[3])
            
            # Use cigar string to adjust for potential soft clipping
            start_pos = soft_checker(line[5])
            
            ## Check bitwise flag for strandedness ##
            bit_flag = bit_checker(int(line[1]))
                
            ## Save chromosome to current (initial) ##
            chromosome = line[2]
            
            ## Save UMI to current (reference) ##
            current_UMI = line[0].split(':')[7]
        	
        	
        	########### Iterate through lines in SAM file to filter out duplicates ############
            if line_cnt > 1:
            
            	# Assumes next line is not a duplicate #
                de_duped = False
                
                ## establishing new variables while calling functions made above ##
                pos = int(line[3])
                next_pos = soft_checker(line[5])
                next_strand = bit_checker(int(line[1]))
                next_chrom = line[2]
            	next_UMI = line[0].split(':')[7]
        
                
                ## Check if UMI is present ##
                if umi_checker(next_UMI) == True:
                	continue
                
                else:
                	if umi_checker(next_UMI) != True
                		de_duped = True
                    	duplicates_cnt +=1
                
                ## Take next line and check parameters and if all 'True', return as duplicate ##
                if chromosome == next_chrom and start_position == next_pos and bit_flag == next_strand and current_UMI == next_UMI:
                    de_duped = True
                    duplicates_cnt += 1
                
                ## if not duplicates found, write out line to file ##
                if de_duped == False:
                    duplicated.write(read_line)
        
                
            ## Reference the first line to see if there is a a valid UMI ##
            if line_cnt == 1:
                if umi_checker(UMI) == True:
                    duplicated.write(read_line)
        	


print("Iteration_complete.")
print("Removed" + str(duplicates_cnt) + " duplications (with invalid UMIs")        	
        	
dulicated.close()       	
        	
        	
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        	
        	
        	
        	
        	
        	

#################################################### Argparse setup ###########################################################

import argparse
import re


parser = argparse.ArgumentParser(description="Iterates through sorted SAM file and removes PCR duplicates from uniq_mapped reads (single & paired end). Quality score per base was used to separate output file while UMIs vs known UMIs were also used to separate output file. Randomers were ignored.")

parser.add_argument('-f', '--file_name', help= "specifying file path for SAM file and read in, r", required=True, type=str)
parser.add_argument('-u', '--UMI', help= "Optional flag designating file containing accepted UMIs from known UMIs list",  required=False, type=str, default='')
parser.add_argument('-p', '--paired_end', help=" File is paired end (not single_end)",required=False, type=bool, default=False)
#parser.add_argument('-h', '--help', help= "prints informative help message about script/output", required=False, type=str, default='')


args = parser.parse_args()


f = args.file_name
u = args.UMI
paired_end = args.paired_end


## When known_UMIs is added ##

if args.UMI:
    known_UMIs = []
    
    with open(args.UMI, 'r+') as umi_document:
        ## for each line in known umi_document ##
        for line in umi_document:
            ## append each UMI to empty list "known_UMIs" ##
            known_UMIs.append(line.strip())

## When known_UMIs is not given ##
else:
    list_not_given = True




## If paired_end option is exercised, exit run with error ##
if args.paired_end:
    raise NameError('Support for paired_end reads is not valid @ this time')

################################# logical grouping of data and functions ######################################

#class sort_sam_file():
    #'''Using class & methods to iterate and filter through SAM file. 
    #Process includes:
        
        #- monitors & adjusts alignment start postion (via soft clipping amount ~ CIGAR string "S").
        #- checks quality score (from conversion of phred score - Using ASCII table ~ assuming ASCII being 33).
        #- checks for strandeness and unique mapping.
        #- checks chromosome positioning
        #- grabs molecular id (umi barcode).
    
    
    #SAM file have the following properties:
        
        #Attributes:
            #line: looking at current_line (single) from sam file (string)
            #start_position: left most mapping postion on said line ~ soft clip adjusted (integer)
            #q_score: line quality score (integer)
            #strand (Bit-wise flag): lines strandeness (string)
            #soft_clipped: determines if seq featured soft clipping (integer)
            #chromosome: which chromosome (integer)
    #'''
    
    
    ## using '_init_' to construct the class (initializing atrributes listed above) and 'self' to call the instance of the object within class ##
    #def __init__(self, current_line):
        
        #self.current_line = line()
        #self.start_position = int(line.split('\t')[3])
        #self.q_score = int(line.split('\t')[4])
        #self.bit_flag = int(line.split('\t')[1])
        #self.soft_clipped = line.split('\t')[5]
        #self.chromosome = int(line.split('\t')[2])
        


## function to check for soft_clipping ##    
def soft_checker(start_position, cigar_string):
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


## function to check for UMI ~ Known_UMI ##
def umi_checker(UMI):
    '''cross-reference current_UMI to known_UMI list ~ site: https://github.com/Leslie-C/Deduper/blob/master/STL96.txt'''
        
    if UMI in known_UMIs:
        return True
    
    else:
        # proceeds to ignore and throw this read out of file # 
        return False


## function to check bit_flag for mapped, strandness, single_vs_paired ##
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
    
    ## looking at paired_end reads; if set 'True'
    if paired_end == True:
        if (bit_flag & 40) != 40 and (bit_flag & 80) == 80:
            ## take first segment ##
            paired_read = 1
        
        elif (bit_flag & 40) == 40 and (bit_flag & 80) != 80: 
            ## take last segment ## 
            paired_read = 2
        
        return paired_read
        


####################################################### Setup for iteration & outputfil format ############################################################

## Set counter to 0 in order to tally up all duplicates from file as well as count for un-used UMIs ## 
duplicates_cnt = 0

## invalid_umi counter ##
invalid_umi_cnt = 0

## set dictionary to hold UMI, pos, strand as key, value as read_line ##
sorted_dic = {}

## set chrom_reference for iteration ##
reference_chrom = ''

## Generate output_sam file w.out duplications in said file ## 
duplicated = open(args.file_name+"_deduped",'w')



####################################################### Iterate through SAM file ############################################################
with open(args.file_name, 'r') as sam_file_1:
    # setting line count to 0 #
    line_cnt = 0
    
    ## looking at each line in SAM file (fh) ##
    for read_line in sam_file_1:
        
        ## iterating through SAM file (HEADER) to see if line starts with certain character (@), if so, ignore and move onto next line##
        if read_line.startswith('@'):
            duplicated.write(read_line)
        
        ## If line consists of read info from SAM file ##
        if read_line.startswith('NS5'):
            ## for each line, increment by 1 ##
            line_cnt += 1
        
            ## looking at each line that is post 'header label' and strips by new line and splits by tab ##
            current_line = read_line.split("\t")
        

    ####################### Initiate reference variables for iteration #######################
            ## Unadjusted raw_pos - soft clipped not checked ##
            pos = int(current_line[3])
            
            ## setting cigar var ##
            cigar_str = str(current_line[5])
            
            ## Use cigar string to adjust for potential soft clipping
            start_pos = soft_checker(pos, cigar_str)
            
            ## Check bitwise flag for strandedness ##
            bit_flag = bit_checker(int(current_line[1]), paired_end)
                
            ## Save chromosome to current (initial) ##
            chromosome = current_line[2]
            
            ## Save UMI to current (reference) ##
            current_UMI = current_line[0].split(':')[7]
            
            ## Check if reference chrom matches current_line ##
            if umi_checker(current_UMI) == True:
                
                ## if chromosome does not match previous chromosome
                if chromosome != reference_chrom:
                ## when reaching new chrom, keep valid lines and output to new file ##
                    for value in sorted_dic.values():
                        duplicated.write(value)
                
                    ## setting new chrom_value ##
                    reference_chrom = chromosome
            
                    ## clear current dic value ##
                    sorted_dic = sorted_dic.clear() 
                    
                    ## reset established dic ##
                    sorted_dic = {}
            
            
                ## setting up input for key position in sorted_dic given UMI being valid and also being on same chrom ##
                Key_input = (current_UMI, start_pos, bit_flag)
                
                ## if dic == {}, add current_line to file (new chrm was hit ##
                if any(sorted_dic) != True: 
                    ## taking that line, given chrom and adding current_line as value ##
                    sorted_dic[Key_input] = read_line 
                    
                ## if dic is not empty ##
                if sorted_dic != {}:
                        
                    ## if the key in dic is pressent, increment by 1 ##
                    if Key_input in sorted_dic: 
                        duplicates_cnt +=1 
                    
                    ## if not present, add read_line to value in sorted_dic ##    
                    else:
                        sorted_dic[Key_input] = read_line
            
            ## If UMI match is not found, dup_cnt += 1 ##
            else:
                if umi_checker(current_UMI) == False:
                    invalid_umi_cnt += 1
                
    ## Given the chromosome being constant through out SAM file ##
    if sorted_dic != {}: 
        for value in sorted_dic.values():
            duplicated.write(value)   



################################### print statements once iteration through SAM file is complete ###############

print("Iteration_complete.")
print("reads_removed" + str(duplicates_cnt) + " duplications (with invalid UMIs") 

duplicated.close()       
            
            


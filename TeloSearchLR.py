# TeloSearchLR
# (TideHunter_tabulated_position_parser_v11o_occupancy_mode.py for distribution)
# uses TideHunter >=v1.5.4 to find tandemly repetitive sequences on reads
# then graphs the occupancies of the different repeats at the terminal n nucleotides of all reads. 

from itertools import chain
from Bio import SeqIO
from PIL import Image
from contextlib import ExitStack
import svgutils.transform as sg
import getopt
import sys
import plotly.graph_objects as go
import os
import math

version_number = "1.0"
TeloSearchLR_output_prefix = "TeloSearchLR_v1"
help_text = """TeloSearchLR: TELOomeric repeat motif SEARCH using Long Reads
                
Version: {0}   Contact: Dr. George Chung (gc95@nyu.edu)
                
Usage:   python {1} [options]
                
Options:
  Run modes:
    (default)                      "occupancy mode", repeat motifs ranked by occupancy
    -e --exhaustive                enable "exhaustive mode", motifs ranked by period AND occupancy
    -s --single_motif       STR    enable "single-motif mode", specify the motif whose occupancy is to be plotted
  Required for all modes:
    -f --fasta              STR    long-read sequencing library file name
    -n --num_of_nucleotides INT    number of nucleotides to plot motif occupancy
  Required for occupancy and exhaustive modes:
    -k --k_value            INT    shortest repeat period (>0) to consider by TideHunter
    -K --K_value            INT    longest repeat period (>=k) to consider by TideHunter
    -m --mth_pattern        INT    most frequent motif to plot (>0)
    -M --Mth_pattern        INT    least frequent motif to plot (>m)
  Required for single-motif mode:
    -T --TideHunter         STR    a TideHunter (>=v1.5.4) tabular output
  Other options:
    -t --terminal           INT    terminal number of nucleotides to consider for ranking motif occupancy [1000]
    -c --cores              INT    number of threads to use for TideHunter [4]
    -p --path               STR    path of TideHunter (if not already in $PATH)
    -v --version                   display the version number and quit
    -h --help                      display this help message and quit\n\n"""

# grab the head and tail ends of nucleotides of every read in a library
def head_and_tail_nucleotides(input_filename, n_nucleotides, output_filename):
    with open(output_filename, "w") as output_fh:
        for current_seq in SeqIO.parse(input_filename, "fasta"):
            if len(current_seq) < 2*n_nucleotides:
                pass
            else:
                #first capture the first n nucleotides, then rename the Seq object id
                current_seq_first_n = current_seq[0:n_nucleotides]
                current_seq_first_n.id = "head_" + str(n_nucleotides) + "_bps_of_" + current_seq_first_n.id
                SeqIO.write([current_seq_first_n], output_fh, "fasta")
                
                #capture the last n nucleotides, then rename the Seq object id
                current_seq_last_n = current_seq[len(current_seq)-n_nucleotides:len(current_seq)]
                current_seq_last_n.id = "tail_" + str(n_nucleotides) + "_bps_of_" + current_seq_last_n.id
                SeqIO.write([current_seq_last_n], output_fh, "fasta")

# extension remover
# removes the last string after the last period of a filename
def filename_extension_remover(filename):
    filename_parts = filename.split(".")
    filename_parts.pop(-1)
    
    return ".".join(filename_parts)

# file path remover
# removes the /path/to/fasta/ if this is supplied in the -f fasta_file.fasta
def filepath_remover(filename):
    filename_parts = filename.split("/")
    
    return filename_parts[-1]

#Define cyclical_permute: a function to circularly permute a string
def cyclical_permute(arg):
	sequence = str(arg)
	return sequence[-1]+sequence[0:len(sequence)-1]

#Define cyclical_permute_matrix function: turns sequence into a list of itself and cyclical permutations
def cyclical_permute_matrix(sequence):
    sequence = sequence.upper()
    return_matrix = [sequence]
    next_cyclical_permutation = cyclical_permute(sequence)
    while next_cyclical_permutation != sequence:
        return_matrix.append(next_cyclical_permutation)
        next_cyclical_permutation = cyclical_permute(return_matrix[len(return_matrix)-1])

    return return_matrix


# Define an index formatter to add 0s to the beginning of indices
# We do not expect to graph more than 999 occupancy patterns
def index_formatter(index):
    index_string =""
    if index < 10:
        index_string="00"+str(index)
    elif index >=10 and index < 100:
        index_string="0"+str(index)
    else:
        index_string=str(index)
    return index_string

#Define reverse_complement: returns a string of the reverse complement of the argument
def reverse_complement(arg):
    revcomp_dictionary = {"A": "T", \
                         "a": "T", \
                         "C": "G", \
                         "c": "G", \
                         "G": "C", \
                         "g": "C", \
                         "T": "A", \
                         "t": "A"}
    sequence = str(arg)
    reverse_complement = ""
    for letter in sequence[::-1]:
        reverse_complement = reverse_complement + revcomp_dictionary[letter]
    
    return reverse_complement

# define head_mid_tail_divider
# Takes an input FASTA file of Nanopore reads (fasta_filename), count the number of nucleotides in the first n_nucleotides
# the last n_nucleotides, and the number of middle nucleotides
def head_mid_tail_divider(fasta_filename, n_nucleotides, head_mid_tail_divider_output_file):
    head_nuc_cumulative = 0
    tail_nuc_cumulative = 0
    midd_nuc_cumulative = 0
    
    num_of_reads_total = 0
    num_of_reads_rejected = 0
    
    for current_seq in SeqIO.parse(fasta_filename, "fasta"):
        if len(current_seq) < 2*n_nucleotides:
            num_of_reads_rejected = num_of_reads_rejected + 1
        else:
            head_nuc_cumulative = head_nuc_cumulative + n_nucleotides
            tail_nuc_cumulative = tail_nuc_cumulative + n_nucleotides
            midd_nuc_cumulative = midd_nuc_cumulative+ len(current_seq) - 2*n_nucleotides
        num_of_reads_total = num_of_reads_total + 1
        
        #string_to_write = "Parsing through " + fasta_filename + ", read #" + str(num_of_reads_total)
        #print(string_to_write, end = "\r")
    
    return_array = [head_nuc_cumulative, midd_nuc_cumulative, tail_nuc_cumulative, num_of_reads_total, num_of_reads_rejected]
    
    with open(head_mid_tail_divider_output_file, "w") as output_fh:
        for value in return_array:
            output_fh.write(str(value) + "\n")
    
    return(return_array)

def head_mid_tail_divider_reader(head_mid_tail_divider_output_file):
    return_array = []
    with open(head_mid_tail_divider_output_file, "r") as output_fh:
        line = output_fh.readline()
        while line:
            return_array.append(int(line))
            line = output_fh.readline()
    
    return(return_array)


# Define TideHunter_tabulated_parse(input_filename, unranked_output_filename, ranked_output_filename)
# Takes the formatted TideHunter tabulated output (input_filename)
# And output 2 files:
#      unranked_output_filename: repeat patterns with their occupancies in the defined space
#      ranked_output_filename: repeat patterns (with their cyclic permutations, and their reverse complements) ranked by total occupancy
def TideHunter_tabulated_parse(input_filename, unranked_output_filename, ranked_output_filename):
    # initialize occupancy_dictionary
    # contains occupancies as counted by TideHunter
    # in the dictionary, key : value pairs are
    #                    repeat : occupancy
    count_dictionary = {}

    # Read the TideHunter output line by line
    with open(input_filename, "r") as filehandle:
        line = filehandle.readline()
        
        while line:
            line_items = line.split("\t")
            
            repeat_pattern = line_items[10].rstrip()
            
            # if the repeat_pattern is already in the dictionary, just add the new TideHunter occupancy value to the current occupancy value
            if repeat_pattern in count_dictionary:
                count_dictionary[repeat_pattern] = count_dictionary[repeat_pattern] + int(line_items[5]) - int(line_items[4]) + 1
            # else, if the repeat_pattern is not already in the dictionary, grab the new TideHunter occupancy value
            else:
                count_dictionary[repeat_pattern] = int(line_items[5]) - int(line_items[4]) + 1
            
            # keep reading through the TideHunter file
            line = filehandle.readline()
    
    '''# Check for / in folder prefix
    if folder_prefix[-1] == "/":
        pass
    else:
        folder_prefix = folder_prefix + "/"
    
    unranked_output_filename = folder_prefix + unranked_output_filename'''

    with open(unranked_output_filename, "w") as filehandle:
        line_to_write = "repeat_pattern\toccupancy\n"
        filehandle.write(line_to_write)
        
        # Now write to file all the repeats and their occupancies
        for key, value in count_dictionary.items():
            line_to_write = key + "\t" + str(value)+ "\n"
            filehandle.write(line_to_write)
    
    # Now, try to consolidate the counts
    # all the cyclical permutations of the same repeat, plus the reverse complements, should be counted as a group
    # ie. ACTG = CTGA = TGAC = GACT (cyclical permutations) = CAGT = AGTC = GTCA = TCAG (reverse complements)
    # data structure: nested list
    #   [ [ pattern_1, occupancy_1], [pattern_2, occupancy_2]...]
    
    consolidated_repeat_counts = []
    
    # step through the count_dictionary and consolidate the counts
    
    while len(count_dictionary) > 0:
        repeat_pattern = list(count_dictionary.keys())[0]
        
        # collect all possible cyclical permutations and reverse complements
        all_cyclical_perms_and_reverse_complements = cyclical_permute_matrix(repeat_pattern) + cyclical_permute_matrix(reverse_complement(repeat_pattern))
        
        # collect and add all occupancy values of each
        total_occupancy = 0
        for pattern in all_cyclical_perms_and_reverse_complements:
            # if the pattern (cyclical permutation or reverse complement) is in the dictionary
            # update the consolidated count
            # pop key (with pattern) from dictionary
            if pattern in count_dictionary:
                total_occupancy = total_occupancy + count_dictionary[pattern]
                count_dictionary.pop(pattern)
            else:
                pass
        
        # store the repeat_pattern and its consolidated counts
        consolidated_repeat_counts.append([repeat_pattern, total_occupancy])
    
    
    #Sort the consolidated_repeat_counts, by the second element (occupancy) of each [repeat, occupancy] pair, descending order
    consolidated_repeat_counts.sort(key = lambda x: x[1], reverse = True)
    
    
    '''ranked_output_filename = folder_prefix + ranked_output_filename
    '''
    
        
    with open(ranked_output_filename, "w") as filehandle:
        # write the consolidated_repeat_counts into a human-readable format
        # with the following columns
        # rank: top ranked occupancies
        # sum_occupancy: sum of occupancies of repeat, cyclical perms, and reverse complements
        # repeat_pattern
        # revcomp_pattern
        # period: period of the repeat in bps
        
        line_to_write = "rank\tsum_occupancy\trepeat_pattern\trevcomp_pattern\tperiod\n"
        filehandle.write(line_to_write)
        
        for index, repeat_occupancy_pair in enumerate(consolidated_repeat_counts):
            line_to_write = str(index+1) + "\t" + str(repeat_occupancy_pair[1]) + "\t" + repeat_occupancy_pair[0] + "\t" + reverse_complement(repeat_occupancy_pair[0]) + "\t" + str(len(repeat_occupancy_pair[0])) + "\n"
            filehandle.write(line_to_write)
    
    return(consolidated_repeat_counts)

# This function reads the file generated by TideHunter_tabulated_parse
def read_parsed_condensed_TideHunter_output(parsed_condensed_TideHunter_output):
    TideHunter_tabulated_parse_matrix = []
    with open(parsed_condensed_TideHunter_output, "r") as parsed_condensed_TideHunter_fh:
        # skip the header line ("rank\tsum_occupancy\trepeat_pattern\trevcomp_pattern\tperiod\n")
        parsed_condensed_TideHunter_fh.readline()
        
        line = parsed_condensed_TideHunter_fh.readline()
        
        
        while line:
            line_items = line.split("\t")
            # capture only the repeat_pattern (line_items[2]) and occupancy (line_items[1])
            TideHunter_tabulated_parse_matrix.append([line_items[2], int(line_items[1])])
            line = parsed_condensed_TideHunter_fh.readline()

    return TideHunter_tabulated_parse_matrix

# split_TideHunter_tabular_output_by_repeat_period(TideHunter_tabular_output_filename, smallest_period, largest_period, split_file_prefix)
# Splits TideHunter output into the -8mer -9mer -10mer etc. files
def split_TideHunter_tabular_output_by_repeat_period(TideHunter_tabular_output_filename, smallest_period, largest_period, split_file_prefix):
    # make a list of filenames
    filename_list = []
    for period in range(smallest_period, largest_period + 1):
        filename_list.append(split_file_prefix + "." + str(period) + "mers.txt")
    
    # open many files at once to be ready to be written
    with ExitStack() as stack:
        files = [stack.enter_context(open(filename, "w")) for filename in filename_list]
        #print(files)
        with open(TideHunter_tabular_output_filename, "r") as TideHunter_tabular_output_fh:
            line = TideHunter_tabular_output_fh.readline()
            
            while line:
                # The 10th item (0-based) of the line is the repeat pattern
                # Thus taking the length of this string would give the period length
                period = len(line.rstrip().split("\t")[10])
                files[period-smallest_period].write(line)
                
                line = line = TideHunter_tabular_output_fh.readline()
    
    #line = "rm " + TideHunter_tabular_output_filename
    #os.system(line)
    


#grapher function
#now returns maximum value on the y-axis (will be important for normalization)
# output files:
# filenames should look like
#   001.repeatPattern.occupancy.head.scatter.txt
#   001.repeatPattern.occupancy.tail.scatter.txt
#   001.repeatPattern.occupancy.all.bars.txt
def head_tail_normalized_grapher(TideHunter_output, num_nucleotides_head_tail, repeat_pattern,
                      head_mid_tail_divider_output, folder_prefix, output_filename_root):
    #Construct a list of repeat_pattern and its circular permutations
    repeat_pattern = repeat_pattern.upper()
    repeat_pattern_matrix = cyclical_permute_matrix(repeat_pattern)
    
    #Define the reverse complement
    revcomp_pattern = reverse_complement(repeat_pattern)
    revcomp_pattern_matrix = []
    
    #construct the revcomp_pattern_matrix only if revcomp_pattern is not already in repeat_pattern_matrix
    if revcomp_pattern in repeat_pattern_matrix:
        pass
    else:
        revcomp_pattern_matrix = cyclical_permute_matrix(revcomp_pattern)
    
    #define and fill the count matrices and fill with num_nucleotides_head_tail x 0s 
    repeat_headcount_matrix = [0]*num_nucleotides_head_tail
    repeat_tailcount_matrix = [0]*num_nucleotides_head_tail
    
    #set the middle count for repeat at 0
    repeat_middcount = 0
    
    #define and fill the count matrices for the reverse complement
    revcomp_headcount_matrix = [0]*num_nucleotides_head_tail
    revcomp_tailcount_matrix = [0]*num_nucleotides_head_tail
    
    #set the middle count for reverse complement at 0
    revcomp_middcount = 0
    
    #num_of_entry = 1

    #print()
    
    #Now open the TideHunter_output to read
    with open(TideHunter_output, "r") as filehandle:
        line = filehandle.readline()
        while line:
            
            #string_to_write = "Examining TideHunter output entry #" + str(num_of_entry) + " for " + repeat_pattern
            #print(string_to_write, end = "\r")
            
            #parse the line read into a list
            line_items = line.split("\t")
            
            #We're interested in
            #line_item[3]: the read length
            #line_item[4]: the start position of a repeat (1-based)
            #line_item[5]: the end position (1-based)
            #line_item[10].rstrip(): the repeat sequence
            read_length = int(line_items[3])
            start_position = int(line_items[4]) - 1
            end_position = int(line_items[5]) - 1
            pattern = line_items[10].rstrip()
            
            #for every line in the TideHunter output, check to see if it matches the repeat pattern(s) in the repeat_pattern_matrix
            #first check to see if the read is long enough (ie. 2x the num_nucleotides_head_tail)
            if read_length < 2*num_nucleotides_head_tail:
                pass

            #else if the read is long enough for analysis
            else:
                if pattern in repeat_pattern_matrix:
                    position_index = start_position
                    
                    while position_index <= end_position:
                        if position_index < num_nucleotides_head_tail:
                            repeat_headcount_matrix[position_index] += 1
                        elif position_index < read_length - num_nucleotides_head_tail:
                            repeat_middcount += 1
                        else:
                            repeat_tailcount_matrix[position_index-read_length+num_nucleotides_head_tail] += 1
                        position_index += 1
                    #now finished filling the repeat_headcount_matrix, repeat_tailcount_matrix, and tallied repeat_middcount
                
                elif pattern in revcomp_pattern_matrix:
                    position_index = start_position
                    
                    while position_index <= end_position:
                        if position_index < num_nucleotides_head_tail:
                            revcomp_headcount_matrix[position_index] += 1
                        elif position_index < read_length - num_nucleotides_head_tail:
                            revcomp_middcount +=1
                        else:
                            revcomp_tailcount_matrix[position_index-read_length+num_nucleotides_head_tail] += 1
                        position_index += 1
                    #now finished filling the revcomp matrices and values
                else:
                    pass
            
            #num_of_entry += 1
            line = filehandle.readline()
    
    #change the revcomp_pattern to "(self-comp/cycl)" if there were no entries in the matrix
    if len(revcomp_pattern_matrix)==0:
        revcomp_pattern = "(self-comp/cycl)"
    else:
        pass
    
    
    #############################################################################################
    # Write to file a summary of repeat/revcomp occupancy at every position at the head         #
    #############################################################################################
    
    if folder_prefix[-1] != "/":
        folder_prefix = folder_prefix + "/"
    else:
        pass
    
    filename = folder_prefix + output_filename_root + ".head.scatter.txt"
    
    # also get the number of reads being considered
    num_of_reads = head_mid_tail_divider_output[3] - head_mid_tail_divider_output[4]
    
    with open(filename, "a") as filehandle:
        string_to_write = "coordinate\tmotif\tcount\toccupancy\n"
        filehandle.write(string_to_write)
        
        # do the repeats first
        for index in range(num_nucleotides_head_tail):
            string_to_write = str(index+1) + "\t" + repeat_pattern + "\t" + str(repeat_headcount_matrix[index]) + "\t"
            filehandle.write(string_to_write)
            
            # change the count to an occupancy
            repeat_headcount_matrix[index] = repeat_headcount_matrix[index]/num_of_reads
            
            # now write the occupancy
            filehandle.write(str(repeat_headcount_matrix[index]) + "\n")
        
        # do the revcomp next
        for index in range(num_nucleotides_head_tail):
            string_to_write = str(index+1) + "\t" + revcomp_pattern + "\t" + str(revcomp_headcount_matrix[index]) + "\t"
            filehandle.write(string_to_write)
            
            # change the count to occupancy
            revcomp_headcount_matrix[index] = revcomp_headcount_matrix[index]/num_of_reads
            
            # now write the occupancy
            filehandle.write(str(revcomp_headcount_matrix[index]) + "\n")

    #############################################################################################
    # Write to file a summary of repeat/revcomp occupancy at every position at the tail         #
    #############################################################################################            
    
    # Wrie to file a summary of repeat/revcomp occupancy at every position at the tail
    
    filename = folder_prefix + output_filename_root + ".tail.scatter.txt"
    
    with open(filename, "a") as filehandle:
        string_to_write = "coordinate\tmotif\tcount\toccupancy\n"
        filehandle.write(string_to_write)
        
        # do the repeats first
        for index in range(num_nucleotides_head_tail):
            string_to_write = str(num_nucleotides_head_tail - index) + "\t" + repeat_pattern + "\t" + str(repeat_tailcount_matrix[index]) + "\t"
            filehandle.write(string_to_write)
            
            # change the count to an occupancy
            repeat_tailcount_matrix[index] = repeat_tailcount_matrix[index]/num_of_reads
            
            # now write the occupancy
            filehandle.write(str(repeat_tailcount_matrix[index]) + "\n")
        
        # do the revcomp next
        for index in range(num_nucleotides_head_tail):
            string_to_write = str(num_nucleotides_head_tail - index) + "\t" + revcomp_pattern + "\t" + str(revcomp_tailcount_matrix[index]) + "\t"
            filehandle.write(string_to_write)
            
            # change the count to occupancy
            revcomp_tailcount_matrix[index] = revcomp_tailcount_matrix[index]/num_of_reads
            
            # now write the occupancy
            filehandle.write(str(revcomp_tailcount_matrix[index]) + "\n")

    #############################################################################################
    # Write to file a summary of repeat occupancy (all regions on the reads) for the bar plot   #
    #############################################################################################  

    filename = folder_prefix + output_filename_root + ".all.bar.txt"
    
    with open(filename, "a") as filehandle:
        filehandle.write("region\tmotif\toccupancy\n")
        
        # repeat_pattern
        string_to_write = "first " + str(num_nucleotides_head_tail) + "bps\t" + repeat_pattern + "\t" \
                            + str(sum(repeat_headcount_matrix)/num_nucleotides_head_tail) + "\n"
        string_to_write = string_to_write + "middle\t" + repeat_pattern + "\t" \
                            + str(repeat_middcount/head_mid_tail_divider_output[1]) + "\n"
        string_to_write = string_to_write + "last " + str(num_nucleotides_head_tail) + "bps\t" + repeat_pattern + "\t" \
                            + str(sum(repeat_tailcount_matrix)/num_nucleotides_head_tail) + "\n"
        
        # revcomp_pattern
        string_to_write = string_to_write + "first " + str(num_nucleotides_head_tail) + "bps\t" + revcomp_pattern + "\t" \
                            + str(sum(revcomp_headcount_matrix)/num_nucleotides_head_tail) + "\n"
        string_to_write = string_to_write + "middle\t" + revcomp_pattern + "\t" \
                            + str(revcomp_middcount/head_mid_tail_divider_output[1]) + "\n"
        string_to_write = string_to_write + "last " + str(num_nucleotides_head_tail) + "bps\t" + revcomp_pattern + "\t" \
                            + str(sum(revcomp_tailcount_matrix)/num_nucleotides_head_tail) + "\n"
        
        filehandle.write(string_to_write)
    
    #############################################################################################
    # figure out the scale for the y-axis taking the maximum of y-values
    # This is 2% more than the maximum of head_count_matrix, revcomp_head_count_matrix,
    # tail_count_matrix, revcomp_tail_count_matrix, repeat_middcount and revcomp_middcount
    #############################################################################################
    
    y_axis_max = 1.02*max(max(repeat_headcount_matrix), max(revcomp_headcount_matrix), \
                            max(repeat_tailcount_matrix), max(revcomp_tailcount_matrix), \
                            repeat_middcount/head_mid_tail_divider_output[1], \
                            revcomp_middcount/head_mid_tail_divider_output[1])
               
    #Now graph the first num_nucleotides_head_tail
    figure_x_axis_index = list(range(1, num_nucleotides_head_tail+1))
    
    
    repeat_name_for_graph = ""
    revcomp_name_for_graph = ""
    
    filename_index = output_filename_root.split(".")[0]
    if len(repeat_pattern) > 37:
        repeat_name_for_graph = "repeat" + filename_index
    else:
        repeat_name_for_graph = repeat_pattern
        
    if len(revcomp_pattern) > 37:
        revcomp_name_for_graph = "revcomp" + filename_index
    else:
        revcomp_name_for_graph = revcomp_pattern
    
    fig_head = go.Figure()
    fig_head.add_trace(go.Scatter(x=figure_x_axis_index, y=repeat_headcount_matrix, mode='lines',
                        name=repeat_name_for_graph))
    fig_head.add_trace(go.Scatter(x=figure_x_axis_index, y=revcomp_headcount_matrix, mode='lines',
                        name=revcomp_name_for_graph))
  
    
    fig_head.update_xaxes(range=[0, len(figure_x_axis_index)])
    fig_head.update_yaxes(range=[0, y_axis_max])
    fig_head.update_layout(
        font_family="Arial",
        autosize=False,
        width=400,
        height=300,
        margin=dict(
            l=70,
            r=50,
            b=50,
            t=50,
            pad=4
        ),
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.05,
            xanchor="right",
            x=1,
            font=dict(
                family="Courier",
                size=12,
                color="black"
            )
        ),
        yaxis = dict(
            showexponent = 'all',
            exponentformat = 'e'
        )
    )
    filename = folder_prefix + output_filename_root + ".head.scatter.svg"
    fig_head.write_image(filename)
    filename = folder_prefix + output_filename_root + ".head.scatter.png"
    fig_head.write_image(filename)
    
    
    #Now graph the last num_nucleotides_head_tail
    figure_x_axis_index.reverse()
    
    
    fig_tail = go.Figure()
    fig_tail.add_trace(go.Scatter(x=figure_x_axis_index, y=repeat_tailcount_matrix, mode='lines',
                        name=repeat_name_for_graph))
    fig_tail.add_trace(go.Scatter(x=figure_x_axis_index, y=revcomp_tailcount_matrix, mode='lines',
                        name=revcomp_name_for_graph))

    fig_tail.update_xaxes(range=[len(figure_x_axis_index), 0])
    fig_tail.update_yaxes(range=[0, y_axis_max])
    fig_tail.update_layout(
        font_family="Arial",
        autosize=False,
        width=400,
        height=300,
        margin=dict(
            l=70,
            r=50,
            b=50,
            t=50,
            pad=4
        ),
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.05,
            xanchor="right",
            x=1,
            font=dict(
                family="Courier",
                size=12,
                color="black"
            )
        ),
        yaxis = dict(
            showexponent = 'all',
            exponentformat = 'e'
        )
    )
    
    filename = folder_prefix + output_filename_root + ".tail.scatter.svg"
    fig_tail.write_image(filename)
    filename = folder_prefix + output_filename_root + ".tail.scatter.png"
    fig_tail.write_image(filename)
    
    #now graph the head, mid, tail fraction occupancy of the repeat class
    
    #define x-axis categories
    figure_x_axis_index = ["first " + str(num_nucleotides_head_tail) + "nts", "middle", "last " + str(num_nucleotides_head_tail) + "nts"]
       
    #define y-axis values for repeat
    repeat_y_values = [sum(repeat_headcount_matrix)/num_nucleotides_head_tail, repeat_middcount/head_mid_tail_divider_output[1], sum(repeat_tailcount_matrix)/num_nucleotides_head_tail]
        
    #define y-axis values for the reverse complement
    revcomp_y_values = [sum(revcomp_headcount_matrix)/num_nucleotides_head_tail, revcomp_middcount/head_mid_tail_divider_output[1], sum(revcomp_tailcount_matrix)/num_nucleotides_head_tail]
           
    #construct bar graph
    fig_bar = go.Figure()
    fig_bar.add_trace(go.Bar(name=repeat_name_for_graph, x=figure_x_axis_index, y=repeat_y_values))
    fig_bar.add_trace(go.Bar(name=revcomp_name_for_graph, x=figure_x_axis_index, y=revcomp_y_values))
    fig_bar.update_yaxes(range=[0, y_axis_max])
    
    fig_bar.update_layout(
        barmode='group',
        font_family="Arial",
        autosize=False,
        width=400,
        height=300,
        margin=dict(
            l=70,
            r=50,
            b=50,
            t=50,
            pad=4
        ),
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.05,
            xanchor="right",
            x=1,
            font=dict(
                family="Courier",
                size=12,
                color="black"
            )
        ),
        yaxis = dict(
            showexponent = 'all',
            exponentformat = 'e'
        )
    )
    filename = folder_prefix + output_filename_root + ".all.bar.svg"
    fig_bar.write_image(filename)
    filename = folder_prefix + output_filename_root + ".all.bar.png"
    fig_bar.write_image(filename)
    
    #now, find the largest y-value for the occupancy graphs
    max_yvalue = max(max(repeat_headcount_matrix), max(repeat_tailcount_matrix), max(revcomp_headcount_matrix), max(revcomp_tailcount_matrix))


# im1, im2 and im3 are PIL.Image objects    
def PNG_image_triple_concat_horizontal(im1, im2, im3):
    dst = Image.new('RGB', (im1.width+im2.width+im3.width, max([im1.height, im2.height, im3.height])))
    dst.paste(im1, (0, 0))
    dst.paste(im2, (im1.width, 0))
    dst.paste(im3, (im1.width+im2.width, 0))
    return dst

# im1 and im2 are PIL.Image objects    
def PNG_image_concat_horizontal(im1, im2):
    dst = Image.new('RGB', (im1.width+im2.width, max([im1.height, im2.height])))
    dst.paste(im1, (0, 0))
    dst.paste(im2, (im1.width, 0))
    return dst

# im1, im2 and im3 are PIL.Image objects
def PNG_image_concat_vertical(im1, im2):
    dst = Image.new('RGB', (max([im1.width, im2.width]), im1.height + im2.height))
    dst.paste(im1, (0, 0))
    dst.paste(im2, (0, im1.height))
    return dst

# im1, im2, im3 are SVG util's filehandles??
def SVG_image_triple_concat_horizontal(im1, im2, im3):
    im1_width = int(im1.get_size()[0])
    im1_height = int(im1.get_size()[1])
    
    im2_width = int(im2.get_size()[0])
    im2_height = int(im2.get_size()[1])
    
    im3_width = int(im3.get_size()[0])
    im3_height = int(im3.get_size()[1])
    
    fig = sg.SVGFigure()    
    fig.set_size((str(im1_width+im2_width+im3_width), str(max(im1_height, im2_height, im3_height))))
    
    # get the plot objects
    im1_root = im1.getroot()
    im2_root = im2.getroot()
    im3_root = im3.getroot()
    
    # move im2 and im3
    im2_root.moveto(im1_width, 0)
    im3_root.moveto(im1_width+im2_width, 0)
    
    # append the fig
    fig.append([im1_root, im2_root, im3_root])
    
    
    return fig

# im1 and im2 are SVG util's filehandles??
def SVG_image_concat_horizontal(im1, im2):
    im1_width = int(im1.get_size()[0])
    im1_height = int(im1.get_size()[1])
    
    im2_width = int(im2.get_size()[0])
    im2_height = int(im2.get_size()[1])
    
    fig = sg.SVGFigure()
    fig.set_size((str(im1_width+im2_width), str(max(im1_height,im2_height))))
    
    # get the plot objects
    im1_root = im1.getroot()
    im2_root = im2.getroot()
    
    # move im2 and im3
    im2_root.moveto(im1_width, 0)
    
    # append the fig
    fig.append([im1_root, im2_root])
    
    return fig

# im1 and im2 are SVG util's filehandles??
def SVG_image_concat_vertical(im1, im2):
    im1_width = int(im1.get_size()[0])
    im1_height = int(im1.get_size()[1])
    
    im2_width = int(im2.get_size()[0])
    im2_height = int(im2.get_size()[1])
    
    fig = sg.SVGFigure()
    fig.set_size((str(max(im1_width, im2_width)), str(im1_height+im2_height)))
    
    # get the plot objects
    im1_root = im1.getroot()
    im2_root = im2.getroot()
    
    # move im2 and im3
    im2_root.moveto(0, im1_height)
    
    # append the fig
    fig.append([im1_root, im2_root])
    
    return fig
    
def main(argv):
    # Start of the script
    
    arg_fasta_file = ""
    arg_k = 0
    arg_K = 0
    arg_m = 0
    arg_M = 0
    arg_n = 0
    arg_t = 1000
    arg_c = "4"
    arg_p = ""
    arg_mode = "occupancy"
    arg_s = ""
    arg_T = ""
    arg_help = help_text.format(version_number, argv[0])
    opts, args = getopt.getopt(argv[1:], "hf:k:K:m:M:n:t:c:p:es:T:v", ["help", "fasta=", "k_value=", "K_value=", "mth_pattern=", "Mth_pattern=", "num_of_nucleotides=", "terminal=", "cores=", "path=", "exhaustive", "single_motif=", "TideHunter=", "version"])
    
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print(arg_help)  # print the help message
            sys.exit(2)
        elif opt in ("-v", "--version"):
            print(version_number)
            sys.exit(2)
        elif opt in ("-f", "--fasta"):
            arg_fasta_file = arg
        elif opt in ("-k", "--k_value"):
            arg_k = int(arg)
        elif opt in ("-K", "--K_value"):
            arg_K = int(arg)
        elif opt in ("-m", "--mth_pattern"):
            arg_m = int(arg)
        elif opt in ("-M", "--Mth_pattern"):
            arg_M = int(arg)
        elif opt in ("-n", "--num_of_nucleotides"):
            arg_n = int(arg)
        elif opt in ("-t", "--terminal"):
            arg_t = int(arg)
        elif opt in ("-c", "--cores"):
            arg_c = arg
            print(arg_c)
        elif opt in ("-p", "--path_of_TideHunter"):
            arg_p = arg
        elif opt in ("-e", "--exhaustive"):
            arg_mode = "exhaustive"
        elif opt in ("-s", "--single_motif"):
            arg_mode = "single motif"
            arg_s = arg
        elif opt in ("-T", "--TideHunter"):
            arg_T = arg

    if arg_mode == "occupancy" or arg_mode == "exhaustive":
        if arg_k <= 0 or arg_K < arg_k or arg_m <= 0 or arg_M < arg_m or arg_n == 0 or arg_fasta_file == "":
            print("Missing or incorrect options. See below.\n")
            print(arg_help)
            sys.exit(2)
        else:
            pass
    elif arg_mode == "single motif":
        if arg_fasta_file == "" or arg_T == "" or arg_n == 0:
            print("Missing or incorrect options. See below.\n")
            print(arg_help)
            sys.exit(2)
        else:
            pass
    else:
        pass
            
    
    if len(arg_p) > 0 and arg_p[-1] != "/":
        arg_p = arg_p + "/"
        
    ###########################################################################################################
    # STEP 1:
    # Generate a FASTA file with just the first and the last t bps (OCCUPANCY/EXHAUSTIVE modes) OR
    # Do nothing (SINGLE MOTIF mode)
    ###########################################################################################################
    head_and_tail_fasta_filename = str(arg_t) + "bp_HT_"+ filepath_remover(arg_fasta_file)
    
    if arg_mode == "occupancy" or arg_mode == "exhaustive":
        print("[main] Step 1: partitioning reads into the first and last " + str(arg_t) + " bps.", end="")
        if os.path.exists(head_and_tail_fasta_filename):
            print(" This file already exists. Skipping.")
            pass
        else:
            head_and_tail_nucleotides(arg_fasta_file, arg_t, head_and_tail_fasta_filename)
            print(" Done.")
    elif arg_mode == "single motif":
        print("[main] Step 1: skipping partitioning reads.")
    
    
    
    ###########################################################################################################
    # STEP 2:
    # Run TideHunter (OCCUPANCY/EXHAUSTIVE modes) OR
    # Check the specified TideHunter file exists (SINGLE MOTIF mode)
    ###########################################################################################################
    
    # check that K is at most half of t
    if arg_k == 0:
        arg_k = 4
    else:
        pass
    
    if arg_K == 0:
        arg_K = math.floor(arg_t/2)
    else:
        pass
    TideHunter_original_reads_output = ""
    TideHunter_head_tail_output = ""
    
    command = arg_p + "TideHunter -v"
    #os.system(command)
    TideHunter_version = os.popen(command).read()
    
    if arg_mode == "occupancy" or arg_mode == "exhaustive":
        print("[main] Step 2: running TideHunter " + TideHunter_version.rstrip() + " to identify tandem repeat motifs.")
        print("[main] Step 2a: runnng TideHunter on the original reads.", end="")
        TideHunter_original_reads_output = TeloSearchLR_output_prefix + ".k" + str(arg_k) + ".K" + str(arg_K) + "." + filepath_remover(filename_extension_remover(arg_fasta_file)) + ".TideHunterTable.txt"
        if os.path.exists(TideHunter_original_reads_output):
            print(" TideHunter run already completed. Skipping.")
            pass
        else:
            command_string = arg_p + "TideHunter -p " + str(arg_k) + " -P " + str(arg_K) + " -f 2 -m " +  str(arg_k) + " -t " + arg_c + " " + arg_fasta_file + " > " + TideHunter_original_reads_output
            os.system(command_string)
            print(" Done.")
        
        print("[main] Step 2b: runnng TideHunter on the partitioned reads.", end="")
        TideHunter_head_tail_output = TeloSearchLR_output_prefix + ".k" + str(arg_k) + ".K" + str(arg_K) + "." + filename_extension_remover(head_and_tail_fasta_filename) + ".TideHunterTable.txt"
        if os.path.exists(TideHunter_head_tail_output):
            print(" TideHunter run already completed. Skipping.")
            pass
        else:
            command_string = arg_p + "TideHunter -p " + str(arg_k) + " -P " + str(arg_K) + " -f 2 -m " +  str(arg_k) + " -t " + arg_c + " " +  head_and_tail_fasta_filename + " > " + TideHunter_head_tail_output
            os.system(command_string)
            print(" Done")
    elif arg_mode == "single motif":
        if os.path.exists(arg_T):
            print("[main] Step 2: Found the specified TideHunter output file.")
        else:
            print("[main] Step 2: Specified TideHunter output file not found. Exiting.")
            sys.exit(2)
        
        
    ##############################################################################################################
    # STEP 3:                                                                                                    #
    # Find the number of bps in the beginning, middle and tail regions in all reads.                             #
    # in head_mid_tail_divider_output (list), the list items are:                                                #
    # [head_nuc_cumulative, midd_nuc_cumulative, tail_nuc_cumulative, num_of_reads_total, num_of_reads_rejected] #
    ##############################################################################################################
    
    print("[main] Step 3: finding the number of bps in the head, middle and tail.", end="")
    head_mid_tail_divider_output = []
    head_mid_tail_divider_file = TeloSearchLR_output_prefix + ".n" + str(arg_n) + "." + filepath_remover(filename_extension_remover(arg_fasta_file)) + ".headTailDivider.txt"
    if os.path.exists(head_mid_tail_divider_file):
        head_mid_tail_divider_output = head_mid_tail_divider_reader(head_mid_tail_divider_file)
        print(" Number of bps in the head/mid/tail determined. Skipping.")
    else:
        #print("\n")
        head_mid_tail_divider_output = head_mid_tail_divider(arg_fasta_file, arg_n, head_mid_tail_divider_file)
        print(" Done.")
    
    #PIL_image_all = Image.new('RGB', (0, 0))
    #SVG_image_all = sg.SVGFigure("0", "0")
    #SVG_image_all.set_size(("0", "0"))


    ###########################################################################################################
    # STEP 4:
    # Reading the TideHunter outputs and graphing
    ###########################################################################################################    
    
    print("[main] Step 4: reading TideHunter counts and graphing repeat occupancies.", end="")
    
    ######################################################
    # running in default sorting: "occupancy" mode       #
    # motifs are sorted by occupancy only                #
    ######################################################
    if arg_mode == "occupancy":
        
        folder_prefix = TeloSearchLR_output_prefix + ".occupancyMode.k" + str(arg_k) + ".K" + str(arg_K) + ".t" + str(arg_t) + "n." + str(arg_n) + "." + filepath_remover(filename_extension_remover(arg_fasta_file)) + ".results"
        
        if not os.path.exists(folder_prefix):
            os.mkdir(folder_prefix)
        else:
            pass
            
        folder_prefix = folder_prefix + "/"
            
        # find out the most frequent kmers from the TideHunter_headtail_1000bps file
        ranked_output_filename = TeloSearchLR_output_prefix + ".k" + str(arg_k) + ".K" + str(arg_K) + ".t" + str(arg_t) + "." + filepath_remover(filename_extension_remover(arg_fasta_file)) + ".rankedRepeatsTable.txt"        
        unranked_output_filename =  TeloSearchLR_output_prefix + ".k" + str(arg_k) + ".K" + str(arg_K) + ".t" + str(arg_t) + "." + filepath_remover(filename_extension_remover(arg_fasta_file)) + ".unrankedRepeatsTable.txt"     
        
        if os.path.exists(ranked_output_filename):
            print(" Ranked repeats table found.", end="")
            TideHunter_parsed_matrix = read_parsed_condensed_TideHunter_output(ranked_output_filename)
        else:
            print(" Generating ranked repeats table.", end="")
            TideHunter_parsed_matrix = TideHunter_tabulated_parse(TideHunter_head_tail_output, unranked_output_filename, ranked_output_filename)
        
            
        PIL_image_0 = Image.new('RGB', (0, 0))
        SVG_image_0 = sg.SVGFigure()
        SVG_image_0.set_size(("0", "0"))
            
        filename_index = arg_m
        
        include_repeat_pattern_in_filenames = False
        include_n_in_filenames = False
        
        while filename_index <= arg_M and filename_index <= len(TideHunter_parsed_matrix):
            repeat_pattern = TideHunter_parsed_matrix[filename_index-1][0].upper()
            # filenames should look like
            #   001.repeatPattern.occupancy.head.scatter.txt
            #   001.repeatPattern.occupancy.tail.scatter.txt
            #   001.repeatPattern.occupancy.all.bars.txt
            output_filename_root = index_formatter(filename_index) + ".repeatPattern.occupancy"
            
            # There are only three text file:s the head file, the tail file, and the bar plot file
            output_head = output_filename_root + ".head.scatter.txt"
            output_tail = output_filename_root + ".tail.scatter.txt"
            output_bars = output_filename_root + ".all.bar.txt"
            
            # include comments for easier debugging later on, as well as including critical info with each run
            output_comments = "# genomic read library: " + arg_fasta_file + "\n"
            output_comments = output_comments + "# repeat pattern: " + repeat_pattern + "\n"
            output_comments = output_comments + "# reverse complement: " + reverse_complement(repeat_pattern) + "\n"
            output_comments = output_comments + "# TeloSearchLR version: " + version_number + "\n"
            output_comments = output_comments + "# TeloSearchLR settings: -k " + str(arg_k) + " -K " + str(arg_K) + " -m " + str(arg_m) + " -M " + str(arg_M) + " -n " + str(arg_n) + " -t " + str(arg_t) + " -f " + arg_fasta_file + "\n"
            output_comments = output_comments + "# TideHunter version: 1.5.4\n"
            output_comments = output_comments + "# TideHunter tabular output: ../" + TideHunter_original_reads_output + "\n"
            output_comments = output_comments + "# ranking table: ../" + ranked_output_filename + "\n"
            output_comments = output_comments + "# number of reads considered: " + str(head_mid_tail_divider_output[3]-head_mid_tail_divider_output[4]) + "\n"
            
            with open(folder_prefix + output_head, "w") as output_fh:
                output_fh.write(output_comments)
                output_fh.write("# occupancy = count/(# of reads considered)\n")
            
            with open(folder_prefix + output_tail, "w") as output_fh:
                output_fh.write(output_comments)
                output_fh.write("# occupancy = count/(# of reads considered)\n")
            
            with open(folder_prefix + output_bars, "w") as output_fh:
                output_fh.write(output_comments)
                output_fh.write("# number of bps in the first " + str(arg_n) + " bps: " + str(head_mid_tail_divider_output[0]) + "\n")
                output_fh.write("# number of bps in the middle: " + str(head_mid_tail_divider_output[1]) + "\n")
                output_fh.write("# number of bps in the last " + str(arg_n) + " bps: " + str(head_mid_tail_divider_output[2]) + "\n")
                
            
            ##############
            
            
            head_tail_normalized_grapher(TideHunter_original_reads_output, arg_n, repeat_pattern,
                          head_mid_tail_divider_output, folder_prefix, output_filename_root)
            
            
            filename = folder_prefix + output_filename_root + ".head.scatter.png"
            PNG_image1 = Image.open(filename)
            
            filename = folder_prefix + output_filename_root + ".tail.scatter.png"
            PNG_image2 = Image.open(filename)
            
            filename = folder_prefix + output_filename_root + ".all.bar.png"
            PNG_image3 = Image.open(filename)
            
            PIL_image_0 = PNG_image_concat_vertical(PIL_image_0, PNG_image_triple_concat_horizontal(PNG_image1, PNG_image2, PNG_image3))
                
            filename = folder_prefix + output_filename_root + ".head.scatter.svg"
            SVG_image1 = sg.fromfile(filename)
                
            filename = folder_prefix + output_filename_root + ".tail.scatter.svg"
            SVG_image2 = sg.fromfile(filename)
                
            filename = folder_prefix + output_filename_root + ".all.bar.svg"
            SVG_image3 = sg.fromfile(filename)
                
            SVG_image_0 = SVG_image_concat_vertical(SVG_image_0, SVG_image_triple_concat_horizontal(SVG_image1, SVG_image2, SVG_image3))
                
                
            filename_index = filename_index + 1
                
            
        filename = folder_prefix + "repeatPattern.m" + str(arg_m) + ".M" + str(arg_M) + ".png"
        PIL_image_0.save(filename)
            
        filename = folder_prefix + "repeatPattern.m" + str(arg_m) + ".M" + str(arg_M) + ".svg"
        SVG_image_0.save(filename)
            
        print("\nScript finished.\n\n")

    #######################################################
    # running in "exhaustive" mode                        #
    # motifs are sorted first by period then by occupancy #
    #######################################################
    
    elif arg_mode == "exhaustive":
        # find out if the split TideHunter output files already exist
        # check TideHunter_head_tail_output
        has_all_the_split_files = True
        
        print("[main] Step 2c: splitting partitioned reads' TideHunter output by period.")
        for period in range(arg_k, arg_K + 1):
            has_all_the_split_files = has_all_the_split_files and os.path.exists(filename_extension_remover(TideHunter_head_tail_output)+ "." + str(period) + "mers.txt")
            if has_all_the_split_files == False:
                split_TideHunter_tabular_output_by_repeat_period(TideHunter_head_tail_output, arg_k, arg_K, filename_extension_remover(TideHunter_head_tail_output))
                print(" Splitting the partitioned TideHunter output by repeat motif period.")
                break
            else:
                pass
        
        if has_all_the_split_files:
            print(" TideHunter output for partitioned reads already split by repeat motif period.")
        else:
            pass
            
        # TideHunter_original_reads_output
        
        # initialize the graphics
        PIL_image_all = Image.new('RGB', (0, 0))
        SVG_image_all = sg.SVGFigure("0", "0")
        SVG_image_all.set_size(("0", "0"))
        
        # set up the folder, eg.
        # TeloSearchLR_v1.exhaustiveMode.k4.K20.t1000.m1.M40.n1000.LIBRARY.results
        folder_prefix = TeloSearchLR_output_prefix + ".exhaustiveMode.k" + str(arg_k) + ".K" \
                        + str(arg_K) + ".t" + str(arg_t) + ".m" + str(arg_m) + ".M" + str(arg_M) \
                        + ".n" + str(arg_n) + "." + filepath_remover(filename_extension_remover(arg_fasta_file)) + ".results"

        if not os.path.exists(folder_prefix):
            os.mkdir(folder_prefix)
        else:
            pass
            
        folder_prefix = folder_prefix + "/"
        
        # plotting
        for period in range(arg_k, arg_K+1):
            
            period_subfolder = folder_prefix + str(period) + "mers"
            
            if not os.path.exists(period_subfolder):
                os.mkdir(period_subfolder)
            else:
                pass
                
            period_subfolder = period_subfolder + "/"
            
            # check for the existence of the ranked and unranked repeats tables:
            # eg. TeloSearchLR_v1.k4.K20.t1000.LIBRARY.unrankedRepeatsTable.4mers.txt
            # find out the most frequent kmers from the TideHunter_headtail_1000bps file
            ranked_output_filename = TeloSearchLR_output_prefix + ".k" + str(arg_k) + ".K" + str(arg_K) + ".t" + str(arg_t) \
                                    + "." + filepath_remover(filename_extension_remover(arg_fasta_file)) + ".rankedRepeatsTable." + str(period) + "mers.txt"
            unranked_output_filename = TeloSearchLR_output_prefix + ".k" + str(arg_k) + ".K" + str(arg_K) + ".t" + str(arg_t) \
                                    + "." + filepath_remover(filename_extension_remover(arg_fasta_file)) + ".unrankedRepeatsTable." + str(period) + "mers.txt"
            
            if os.path.exists(ranked_output_filename):
                print(" Ranked repeats table found for k = " + str(period) + ".", end="")
                TideHunter_parsed_matrix = read_parsed_condensed_TideHunter_output(ranked_output_filename)
            else:
                print(" Generating ranked repeats table for k = " + str(period) + ".", end="")
                TideHunter_parsed_matrix = TideHunter_tabulated_parse(filename_extension_remover(TideHunter_head_tail_output)+ "." + str(period) + "mers.txt", unranked_output_filename, ranked_output_filename)
            
            PIL_image_0 = Image.new('RGB', (0, 0))
            SVG_image_0 = sg.SVGFigure()
            SVG_image_0.set_size(("0", "0"))
            
            filename_index = arg_m
            
            #TideHunter_filename = filename_extension_remover(TideHunter_original_reads_output)+ "." + str(period) + "-mers.txt"
            
            print(" Plotting...")
            
            while filename_index <= arg_M and filename_index <= len(TideHunter_parsed_matrix):
                repeat_pattern = TideHunter_parsed_matrix[filename_index-1][0].upper()
                
                # filenames should look like
                #   001.repeatPattern.occupancy.head.scatter.txt
                #   001.repeatPattern.occupancy.tail.scatter.txt
                #   001.repeatPattern.occupancy.all.bars.txt
                output_filename_root = index_formatter(filename_index) + ".repeatPattern.occupancy"
                
                # There are only three text file:s the head file, the tail file, and the bar plot file
                output_head = output_filename_root + ".head.scatter.txt"
                output_tail = output_filename_root + ".tail.scatter.txt"
                output_bars = output_filename_root + ".all.bar.txt"
                
                # include comments for easier debugging later on, as well as including critical info with each run
                output_comments = "# genomic read library: " + arg_fasta_file + "\n"
                output_comments = output_comments + "# repeat pattern: " + repeat_pattern + "\n"
                output_comments = output_comments + "# reverse complement: " + reverse_complement(repeat_pattern) + "\n"
                output_comments = output_comments + "# TeloSearchLR version: " + version_number + "\n"
                output_comments = output_comments + "# TeloSearchLR settings: -k " + str(arg_k) + " -K " + str(arg_K) + " -m " + str(arg_m) + " -M " + str(arg_M) + " -n " + str(arg_n) + " -t " + str(arg_t) + " -f " + arg_fasta_file + "\n"
                output_comments = output_comments + "# TideHunter version: 1.5.4\n"
                output_comments = output_comments + "# TideHunter tabular output: ../" + TideHunter_original_reads_output + "\n"
                output_comments = output_comments + "# ranking table: ../" + ranked_output_filename + "\n"
                output_comments = output_comments + "# number of reads considered: " + str(head_mid_tail_divider_output[3]-head_mid_tail_divider_output[4]) + "\n"
                
                with open(period_subfolder + output_head, "w") as output_fh:
                    output_fh.write(output_comments)
                    output_fh.write("# occupancy = count/(# of reads considered)\n")
                
                with open(period_subfolder + output_tail, "w") as output_fh:
                    output_fh.write(output_comments)
                    output_fh.write("# occupancy = count/(# of reads considered)\n")
                
                with open(period_subfolder + output_bars, "w") as output_fh:
                    output_fh.write(output_comments)
                    output_fh.write("# number of bps in the first " + str(arg_n) + " bps: " + str(head_mid_tail_divider_output[0]) + "\n")
                    output_fh.write("# number of bps in the middle: " + str(head_mid_tail_divider_output[1]) + "\n")
                    output_fh.write("# number of bps in the last " + str(arg_n) + " bps: " + str(head_mid_tail_divider_output[2]) + "\n")
                
                head_tail_normalized_grapher(TideHunter_original_reads_output, arg_n, repeat_pattern,
                          head_mid_tail_divider_output, period_subfolder, output_filename_root)
                
                
                filename = period_subfolder + output_filename_root + ".head.scatter.png"
                PNG_image1 = Image.open(filename)
            
                filename = period_subfolder + output_filename_root + ".tail.scatter.png"
                PNG_image2 = Image.open(filename)
            
                filename = period_subfolder + output_filename_root + ".all.bar.png"
                PNG_image3 = Image.open(filename)
            
                PIL_image_0 = PNG_image_concat_vertical(PIL_image_0, PNG_image_triple_concat_horizontal(PNG_image1, PNG_image2, PNG_image3))
                
                filename = period_subfolder + output_filename_root + ".head.scatter.svg"
                SVG_image1 = sg.fromfile(filename)
                
                filename = period_subfolder + output_filename_root + ".tail.scatter.svg"
                SVG_image2 = sg.fromfile(filename)
                
                filename = period_subfolder + output_filename_root + ".all.bar.svg"
                SVG_image3 = sg.fromfile(filename)
                
                SVG_image_0 = SVG_image_concat_vertical(SVG_image_0, SVG_image_triple_concat_horizontal(SVG_image1, SVG_image2, SVG_image3))
                
                
                filename_index = filename_index + 1
                
            
            filename = period_subfolder + "repeatPattern.m" + str(arg_m) + ".M" + str(arg_M) + ".png"
            PIL_image_0.save(filename)
            PIL_image_all = PNG_image_concat_horizontal(PIL_image_all, PIL_image_0)
            
            
            filename = period_subfolder + "repeatPattern.m" + str(arg_m) + ".M" + str(arg_M) + ".svg"
            SVG_image_0.save(filename)
            SVG_image_all = SVG_image_concat_horizontal(SVG_image_all, SVG_image_0)
            
            #line_to_write = "sed -i 's/\&\#226;\&\#710;\&\#8217;/-/g' " + filename
            #print(line_to_write)
            #os.system(line_to_write)
            
        
        all_image_file_root = folder_prefix + "repeatPattern.m" + str(arg_m) + ".M" + str(arg_M)
                
        PIL_image_all.save(all_image_file_root + ".png")
        SVG_image_all.save(all_image_file_root + ".svg")
            
        print("\nScript finished.\n\n")
        
    elif arg_mode == "single motif":
        folder_prefix = TeloSearchLR_output_prefix + ".singleMotifMode." + filepath_remover(filename_extension_remover(arg_fasta_file)) + ".results"
        if not os.path.exists(folder_prefix):
            os.mkdir(folder_prefix)
        else:
            pass
        
        repeat_pattern = arg_s
        folder_prefix = folder_prefix + "/"
        output_filename_root = repeat_pattern + ".occupancy.n" + str(arg_n)
        
        # filenames should look like
        #   TTAGGC.occupancy.n1000.head.scatter.txt
        #   TTAGGC.occupancy.n1000.occupancy.tail.scatter.txt
        #   TTAGGC.occupancy.n1000.all.bar.txt
        output_filename_root = repeat_pattern + ".occupancy.n" + str(arg_n) 
        
        # There are only three text file:s the head file, the tail file, and the bar plot file
        output_head = output_filename_root + ".head.scatter.txt"
        output_tail = output_filename_root + ".tail.scatter.txt"
        output_bars = output_filename_root + ".all.bar.txt"
            
        # include comments for easier debugging later on, as well as including critical info with each run
        output_comments = "# genomic read library: " + arg_fasta_file + "\n"
        output_comments = output_comments + "# repeat pattern: " + repeat_pattern + "\n"
        output_comments = output_comments + "# reverse complement: " + reverse_complement(repeat_pattern) + "\n"
        output_comments = output_comments + "# TeloSearchLR version: " + version_number + "\n"
        output_comments = output_comments + "# TeloSearchLR settings: -f " + arg_fasta_file + " -n " + str(arg_n) + " -s " + arg_s + " -T " + arg_T + "\n"
        output_comments = output_comments + "# TideHunter tabular output: ../" + arg_T + "\n"
        output_comments = output_comments + "# number of reads considered: " + str(head_mid_tail_divider_output[3]-head_mid_tail_divider_output[4]) + "\n"
        
        with open(folder_prefix + output_head, "w") as output_fh:
            output_fh.write(output_comments)
            output_fh.write("# occupancy = count/(# of reads considered)\n")
            
        with open(folder_prefix + output_tail, "w") as output_fh:
            output_fh.write(output_comments)
            output_fh.write("# occupancy = count/(# of reads considered)\n")
            
        with open(folder_prefix + output_bars, "w") as output_fh:
            output_fh.write(output_comments)
            output_fh.write("# number of bps in the first " + str(arg_n) + " bps: " + str(head_mid_tail_divider_output[0]) + "\n")
            output_fh.write("# number of bps in the middle: " + str(head_mid_tail_divider_output[1]) + "\n")
            output_fh.write("# number of bps in the last " + str(arg_n) + " bps: " + str(head_mid_tail_divider_output[2]) + "\n")
                
            
        ##############
        
        head_tail_normalized_grapher(arg_T, arg_n, repeat_pattern,
                                    head_mid_tail_divider_output, folder_prefix, output_filename_root)
           
        PIL_image_0 = Image.new('RGB', (0, 0))
        SVG_image_0 = sg.SVGFigure()
        SVG_image_0.set_size(("0", "0"))
        
        filename = folder_prefix + output_filename_root + ".head.scatter.png"
        PNG_image1 = Image.open(filename)
            
        filename = folder_prefix + output_filename_root + ".tail.scatter.png"
        PNG_image2 = Image.open(filename)
            
        filename = folder_prefix + output_filename_root + ".all.bar.png"
        PNG_image3 = Image.open(filename)
            
        PIL_image_0 = PNG_image_triple_concat_horizontal(PNG_image1, PNG_image2, PNG_image3)
        filename = folder_prefix + output_filename_root + ".all.collage.png"
        PIL_image_0.save(filename)
        
                
        filename = folder_prefix + output_filename_root + ".head.scatter.svg"
        SVG_image1 = sg.fromfile(filename)
                
        filename = folder_prefix + output_filename_root + ".tail.scatter.svg"
        SVG_image2 = sg.fromfile(filename)
                
        filename = folder_prefix + output_filename_root + ".all.bar.svg"
        SVG_image3 = sg.fromfile(filename)
                
        SVG_image_0 = SVG_image_triple_concat_horizontal(SVG_image1, SVG_image2, SVG_image3)
        filename = folder_prefix + output_filename_root + ".all.collage.svg"
        SVG_image_0.save(filename)
        
        print("\nScript finished.\n\n")
        
    
    
if __name__ == "__main__":
    main(sys.argv)
    




		
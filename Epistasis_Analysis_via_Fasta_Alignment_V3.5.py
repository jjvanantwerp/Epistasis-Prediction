# Epistasis Analysis via Fasta Alignment
# This is a program that looks for potential epistasis in a fasta alignment.
# Written by James VanAntwerp in May 2020
# Written for the Woldring Lab, Michigan State University in East Lansing, Michigan, USA.

import Bio
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from itertools import combinations
from math import log2
import csv
import numpy
import matplotlib.mlab as mlab
import matplotlib.pyplot as pyplot
from statistics import stdev
import os

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# These are the lists of information that we're going to want to have at different points on a global scale.
High_RSE_list =[]
Low_RSE_list =[]
AAs_in_each_position =[]
Significant_no_pairing_list = []
OATP_Sequences_list =[]
hist_plot_list = []
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~Open the control file and read in settings:~~~~~~~~~~
# From revision one to revision two, the ordering of these items has changed, and it no longer takes in upper/lower bound.
# From version two to version 3.5, this now has a defined upper/lower value.
infile_name=""
upper_bound=-0.1
lower_bound=0.1
outfile_name=""
# Open the control file, and take in all the lines as a list of stirngs.
with open('control.txt') as control:
    lines = list(control)
    # Read in over all the lines. We'll use some string handling to stip them down to just the info we want 
    # (no leading/trailing spaces, just numbers, etc), and it will handle 2 or 4 line input files. Users can
    # manually enter their high/low values, otherwise we'll set them later using standard deviation.
    for i in range (len(lines)):
        num = lines[i].find(":")
        # line two: infile - mandatory
        if (i==1):
            infile_name = (lines[i])[num+1:-1].strip()
            print(f"Input taken from {infile_name}")
        # Line three: outfile - mandatory
        if (i==2):
            if ((lines[i])[-1]=='\n'):
                outfile_name = ((lines[i])[num+1:-1].strip())
            else:
                outfile_name = ((lines[i])[num+1:].strip())
            print (f"Outfile: {outfile_name}")
        # Line four - upper bound - optional
        if ((i==3)and(lines[i]!="\n")):
            upper_bound_temp = (lines[i])[num+1:-1].strip()
            upper_bound = float(upper_bound_temp)
            print(f"Upper bound: {upper_bound}")
       # Line five - lower bound - optional
        if ((i==4)and(lines[i]!="\n")) :
            lower_bound_temp = (lines[i])[num+1:-1].strip()
            lower_bound = float(lower_bound_temp)
            print(f"Lower bound: {lower_bound}")  
# We now have the control conditions set up.

#~~~~~~~~~~Open the input file and read in fasta data:~~~~~~~~~~
# We assume here that the infile is a fasta-formated alignment
with open(infile_name) as handle:
    for Seq_Rec_Itr in SeqIO.parse(handle, "fasta"):
        OATP_Sequences_list.append(str(Seq_Rec_Itr.seq))
num_sequences = len(OATP_Sequences_list)
length_sequence = len(OATP_Sequences_list[0])
# We now have a list of all our sequences in the alignment. Each item in OATP_Sequences_list is a string that is a sequence, and each 
# item in that string is an alligned postion, with '-' standing in for gaps.

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def RSE_Function (AA1,AA2,pair):
    # We take in two amino acids, and a pair from our combinatorial funciton that tells us which two posiitons we're examining.
    AA_count_1 =0
    AA_count_2 =0
    EOPF_count = 0
    b1 = False
    b2 = False
    # Set up some variables
    for sequence in range (num_sequences): # For each sequence, at the given positions. We go down the list of seuqences, and look at each collumn.
        if ( ((OATP_Sequences_list[sequence]) [pair[0]])  ==AA1):
            # If AA1 is in position one
            AA_count_1+=1
            b1=True # We did find an amino acid in this sequence
        if ( ((OATP_Sequences_list[sequence]) [pair[1]])  ==AA2):
            # If AA2 is in position two
            AA_count_2+=1
            b2=True # We did find an amino acid in this sequence
        if(b1 and b2):
            # If both AA1 and AA2 occur in this sequence, add to EOPFcount
            EOPF_count+=1
        # reset our booleans 
        b1=False
        b2=False
    # we should now have a count for the number of instances of AA1, AA2, and Both in the sequences list
    freq_1 = float(AA_count_1)/float(num_sequences)
    freq_2 =float(AA_count_2)/float(num_sequences)
    EPOF = float(EOPF_count)/float(num_sequences)
    # If the two occured together in at least one sequence, our math is striagtforward; return the RSE value
    if (EPOF!=0):
        return (EPOF * log2(EPOF/(freq_1*freq_2)))
    # If the two amino acids do NOT occur together, we need to determine if that's significant. Maybe they're single mutants, 
    # or maybe they're large contributors to the site. We'll figure out how often they SHOULD be together with thier individual frequencies.
    # If their lack of congruence is significant, let's say something.
    #elif((freq_1*freq_2)>0.005):
    #    return "Significant Disunion"
    else:
        return "Nothing Useful"  
# This function will tell us the likelyhood of an epistatic relationship between two sites. It's very good at noticing and quantifying positive
# correlations, but has to infer where there are "holes" in the data.

#~~~~~~~~~~Create a list of lists of the amino acids found in each position:~~~~~~~~~~
for position in range (length_sequence):
    AA_in_pos_temp = []
    # First go to a positon, then go down the column and tally unique ammino acids
    for sequence in range (num_sequences):
        if ((OATP_Sequences_list[sequence])[position] not in AA_in_pos_temp):
            AA_in_pos_temp.append((OATP_Sequences_list[sequence])[position])
    AAs_in_each_position.append(AA_in_pos_temp)
# We now have a list the length of our alignment which contains a list of unique amino acids for each positon. This will help ensure we check all 
# possibilites for epistatic relationships without duplicating a call to our RSE function.

#~~~~~~~~~~Taking RSE values~~~~~~~~~~                
# This section is a little complex. If the user provided upper and lower bounds, it's straightfoward (See the else statment below for that code), 
# but otherwise this becomes much harder. If no upper/lower bound was provided, we're going to use tiwce the standard distibution of the RSE value.
# The dificulty arrises because until we actually run the whole analysis, we don't know the RSE values. So instead, we'll go through the fist 15% of
# the sample, storing all the data. Once we are 15% completed, we estimate the standard distribution for the RSE values of the whole data set,
# and sort the data we have been previously storing into high/low lists. We can then continue the last 90% of the run more efficently. It is tempting
# to wait for the whole dataset to determine the standard distribution, but it can be computationally prohibitive to store all 280,000 data points
# that we generate, and then to itterate over that list a second time. 

# If the user has not provided upper/lower bounds
if ((upper_bound == 0) and (lower_bound == 0)):
    break_bool = True
    combinations_data_storage_list =[]
    for pair in list(combinations(range(length_sequence),2)):      
    # Section one - first 15%
        if (pair[0]<(length_sequence*.15)):
            for AA1 in AAs_in_each_position[pair[0]]: # For each unique amino acid in position one
                for AA2 in AAs_in_each_position[pair[1]]: # For each unique amino acid in position two
                    # This will loop thorugh each combination of unique amino acids from position one and position two.
                    value = RSE_Function(AA1,AA2,pair)
                    # If the RSE value is a number, add it to hist_plot list. Store the data for later evaluation of high/low.
                    if (not isinstance(value, str)):
                        hist_plot_list.append(value)
                        combinations_data_storage_list.append([AA1,pair[0],AA2,pair[1],value])
                    # If our RSE function has found a significant disunion, write it down.
                    if (value == "Significant Disunion"):
                        Significant_no_pairing_list.append([AA1,pair[0],AA2,pair[1]])                     
     # Last 85%
        else:
            # Now let's assign high and low RSE values, but only once
            if(break_bool):
                break_bool = False # This ensures we only run this block of code once
                # Find the standard deviation of our hist_plot_list so far, and we'll take two standard deviations (5%) as interesting.
                RSE_stdev = stdev(hist_plot_list)
                lower_bound = RSE_stdev * -2
                upper_bound = RSE_stdev * 2
                # Now let's sort the data we stored from the first 15% of out run.
                for data_point in (combinations_data_storage_list):    
                    if (data_point[4]>lower_bound):
                        High_RSE_list.append(data_point)
                    if (data_point[4]<upper_bound):
                        Low_RSE_list.append(data_point)
            # We can now continue for the last 85% of our pairs.
            for AA1 in AAs_in_each_position[pair[0]]: # For each unique amino acid in position one
                for AA2 in AAs_in_each_position[pair[1]]: # For each unique amino acid in position two
                    # This will loop thorugh each combination of unique amino acids from position one and position two.
                    value = RSE_Function(AA1,AA2,pair)
                    # If the RSE value is a number, sort into high/low, after adding RSE value to hist_plot list.
                    if (not isinstance(value, str)):
                        hist_plot_list.append(value)
                        if (value>lower_bound):
                            High_RSE_list.append([AA1,pair[0],AA2,pair[1],value])
                        if (value<upper_bound):
                            Low_RSE_list.append([AA1,pair[0],AA2,pair[1],value])
                    # If our RSE function has found a significant disunion, write it down.
                    elif (value == "Significant Disunion"):
                        Significant_no_pairing_list.append([AA1,pair[0],AA2,pair[1]])   

# If the user has provided upper/lower bounds
else:
    Num_Comparisons = 0
    for pair in list(combinations(range(length_sequence),2)):
        # This gives us a for loop that will itterate over every possible two-membered combination of position in out OATP_Sequences_list
        # As of 5:30PM on May 11, the for loop correctly outputs ~283,000 non-duplicate pairs of numbers from 0-753
        for AA1 in AAs_in_each_position[pair[0]]: #for each unique amino acid in position one
            for AA2 in AAs_in_each_position[pair[1]]: #for each unique amino acid in position two
                Num_Comparisons +=1
                # This will loop thorugh each combination of unique amino acids from position one and position two.
                value = RSE_Function(AA1,AA2,pair)
                if (not isinstance(value, str)):
                    hist_plot_list.append(value)
                    if (value>lower_bound):
                        High_RSE_list.append([AA1,pair[0],AA2,pair[1],value])
                    if (value<upper_bound):
                        Low_RSE_list.append([AA1,pair[0],AA2,pair[1],value])
                #elif (value == "Significant Disunion"):
                #    Significant_no_pairing_list.append([AA1,pair[0],AA2,pair[1]])
    print("Number of comparisions: %s" % Num_Comparisons)

#~~~~~~~~~~Making a directory for storing our output files~~~~~~~~~~ 
try:
    # I'm not confident about the error-handling of this whole block, so if it fails, we'll have to abort.
    path = os.getcwd()
    adendum = ""
    i=0
    # This loop will ensure that out directory name is valid, and add a number to the end (1,2,3,4,etc.) util we get a valid directory name
    while (os.path.exists(f"./{outfile_name}{adendum}")):
        i+=1
        adendum=str(i)
    dirout_name = outfile_name+adendum
    os.mkdir(dirout_name)
    print(f"Directory for output made: {dirout_name}")
except:
    print("There was an error creating the directory for output files. Please try again.")

#Let's make a few other useful data structures
High_RSE_Heatmap = numpy.zeros((length_sequence,length_sequence))
Low_RSE_Heatmap = numpy.zeros((length_sequence,length_sequence))
Disunion_RSE_Heatmap = numpy.zeros((length_sequence,length_sequence))

#~~~~~~~~~~Writing AAs, position, and HIGH RSE values to CSV files~~~~~~~~~~                   
with open(f"./{dirout_name}/{outfile_name}_High_Values.csv",mode="w+") as csv_out:
    pen=csv.writer(csv_out,delimiter=',')
    # Header
    pen.writerow(["AA1","Position 1","AA2","Position 2","RSE Value"])
    # All the data
    for item in High_RSE_list:
        pen.writerow([item[0],item[1]+1,item[2],item[3]+1,item[4]])
        ((High_RSE_Heatmap[item[1]]) [item[3]]) +=1 # Acess the position in the heatmap represneting this pairing of rows, and add one to the heatmap.
numpy.savetxt(f"./{dirout_name}/{outfile_name}_High_RSE_Heatmap.csv", High_RSE_Heatmap, fmt='%04d', delimiter=",")

#~~~~~~~~~~Writing AAs, position, and LOW RSE values to CSV file~~~~~~~~~~
with open(f"./{dirout_name}/{outfile_name}_Low_Values.csv",mode="w+") as csv_out:
    pen=csv.writer(csv_out,delimiter=',')
    # Header
    pen.writerow(["AA1","Position 1","AA2","Position 2","RSE Value"])
    # All the data
    for item in Low_RSE_list:
        pen.writerow([item[0],item[1]+1,item[2],item[3]+1,item[4]])  
        ((Low_RSE_Heatmap[item[1]]) [item[3]]) +=1 # Acess the position in the heatmap represneting this pairing of rows, and add one to the heatmap.
numpy.savetxt(f"./{dirout_name}/{outfile_name}_Low_RSE_Heatmap.csv", Low_RSE_Heatmap, fmt='%04d', delimiter=",")

#~~~~~~~~~~Writing AAs, position, and pairing DISUNION values to CSV file~~~~~~~~~~
with open(f"./{dirout_name}/{outfile_name}_Disunion.csv",mode="w+") as csv_out:
    pen=csv.writer(csv_out,delimiter=',')
    # Header
    pen.writerow(["AA1","Position 1","AA2","Position 2"])
    # All the data
    for item in Significant_no_pairing_list:
        pen.writerow([item[0],item[1]+1,item[2],item[3]+1]) 
        ((Disunion_RSE_Heatmap[item[1]]) [item[3]]) +=1 # Acess the position in the heatmap represneting this pairing of rows, and add one to the heatmap.
numpy.savetxt(f"./{dirout_name}/{outfile_name}_Disunion_RSE_Heatmap.csv", Disunion_RSE_Heatmap, fmt='%04d', delimiter=",")

#~~~~~~~~~~Make a histogram!~~~~~~~~~~
pyplot.hist(hist_plot_list, bins=200)
pyplot.xlabel("HSE Value")
pyplot.ylabel("Number of pairings")
pyplot.title("Hisotgram of OATP Epistasis Values")
pyplot.savefig(f"./{dirout_name}/{outfile_name}")
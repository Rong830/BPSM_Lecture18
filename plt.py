#!/usr/local/bin/python3

'''
The files required for today's exercise are on the server, in the directory
/localdisk/data/BPSM/Lecture18/
Use the ecoli.txt file and, using the code above as a starting point, make a chart which shows the AT content in a sliding window.
The alignment.txt file contains a multiple sequence alignment, one sequence per line.
 We want to draw a plot of protein conservation based on the alignment similarities, but without using plotcon!

Your programme/script will have to :
read in the data
get each column in turn
measure the conservation at each column
append the measurement to a list
draw the chart
The interesting bit is how to measure the conservation. There are lots of ways, just chose whichever one you want, for example:
how many unique aa residues in this column
what proportion of aa residues are the same
1 minus sum of squares of proportion of each residue
something more complicated
'''

import matplotlib.pyplot as plt
import numpy as np


# Open the file, read it, and remove the newlines to get one long string
# then take the first 100000 characters
ecoli = open("/localdisk/data/BPSM/Lecture18/ecoli.txt").read().replace('\n', '').upper()[0:100000]

# Set sliding window size
window = 1000

# Four lists to hold the numbers of each base
a = []
t = []
g = []
c = []

# iterate over all the starting positions
for start in range(len(ecoli) - window):
# get the current sliding window
    win = ecoli[start:start+window]
# count each of the four bases and append to the list
    a.append(win.count('A') / window)
    t.append(win.count('T') / window)
    g.append(win.count('G') / window)
    c.append(win.count('C') / window)



# Use the ecoli.txt file and, using the code above as a starting point, make a chart which shows the AT content in a sliding window.
# Set sliding window size
window = 1000
# Four lists to hold the numbers of each base
AT_content = []
# iterate over all the starting positions
for start in range(len(ecoli) - window):
# get the current sliding window
    win = ecoli[start:start+window]
# count each of the four bases and append to the list
    AT_content.append(win.count('A') + win.count('T') /win.count('A') + win.count('T') + win.count('G') + win.count('C'))


plt.figure(figsize=(20,10))
plt.title('AT content in the ecoli sequences')
plt.plot(AT_content)
plt.show()



# The alignment.txt file contains a multiple sequence alignment, one sequence per line.
# We want to draw a plot of protein conservation based on the alignment similarities, but without using plotcon!
file = open("/localdisk/data/BPSM/Lecture18/alignment.txt")
align = file.readlines()
align = [i.replace('\n','') for i in align]
similar = []
for i in range(len(align[0])):
    count = 0
    for a in align:
        if a[i] != '-':
            count +=1
    similar.append(count/len(align[0]))    


plt.figure(figsize=(20,10))
plt.title('protein conservation based on the alignment similarities')
plt.plot(similar)
plt.show()


# how many unique aa residues in this column
file = open("/localdisk/data/BPSM/Lecture18/alignment.txt")
align = file.readlines()
align = [i.replace('\n','') for i in align]
unique_aa = []
for i in range(len(align[0])):
    count = []
    for a in align:
        if a[i] != '-':
            count.append(a[i])
    unique_aa.append(len(set(count)))  


plt.figure(figsize=(20,10))
plt.title('unique aa residues')
plt.plot(unique_aa)
plt.show()

# what proportion of aa residues are the same
file = open("/localdisk/data/BPSM/Lecture18/alignment.txt")
align = file.readlines()
align = [i.replace('\n','') for i in align]
aa = []
for i in range(len(align[0])):
    count = []
    for a in align:
        if a[i] != '-':
            count.append(a[i])
    aa.append(len(set(count))/len(count))  


plt.figure(figsize=(20,10))
plt.title('unique aa residues')
plt.plot(aa)
plt.show()

# 1 minus sum of squares of proportion of each residue
align = open("/localdisk/data/BPSM/Lecture18/alignment.txt").read().replace('\n', '').upper()
ss = []
label = []
count = []
for i in set(list(align)):
    count.append(list(align).count(i))
    label.append(i)


mean = np.mean(count)
for i in count:
    ss.append(1-(i-mean)**2)


plt.figure(figsize=(20,10))
plt.xticks=label
plt.title('1 minus sum of squares of proportion of each residue')
plt.plot(ss)
plt.show()




# Make a chart which shows the AT content in a sliding window
# Written by s123456 on 19 Nov 2021
import os
import matplotlib.pyplot as plt
# Open the file, read it, and remove the newlines to get one long string
# then take the first 100000 characters
ecoli = open("/localdisk/data/BPSM/Lecture18/ecoli.txt").read().replace('\n', '').upper()[0:100000]
# Set sliding window size
window = 1000
# A list to hold the numbers
at = []
# Iterate over all the starting positions
for start in range(len(ecoli) - window):
    win = ecoli[start:start+window]
    at.append(  (win.count('A')+win.count('T')) / window)


# Plot the list with appropriate labels
plt.figure(figsize=(20,10))
plt.plot(at, label="AT content",linewidth=3,color="purple")
plt.ylabel('Fraction of bases')
plt.xlabel('Position on genome')
plt.title("EXERCISE 1\nAT content, 1kb windows of the E.coli genome")
plt.legend()
plt.savefig("Chart_Exercise_1.png",transparent=True)
plt.show()


# Read the alignment into a list of strings
aln = open("/localdisk/data/BPSM/Lecture18/alignment.txt")
aligned_seqs = []
counter =0
# check what length they are as we go
for line in aln:
    counter += 1
    print("Sequence",counter,"was",len(line.rstrip("\n")),"long")
    aligned_seqs.append(line.rstrip("\n"))

alignment_length = len(aligned_seqs[0])
uniques_per_column = []

# Look at each column
for column_number in range(alignment_length):
    column = []
    for seq in aligned_seqs:
        aa = seq[column_number]
        if aa != '-':  # ignore gaps
            column.append(seq[column_number])
    uniques = len(set(column))
    uniques_per_column.append(uniques)
    

# Now do a sliding window over the numbers
window = 10
numbers_to_plot = []
for start in range(len(uniques_per_column) - window):
    win = uniques_per_column[start:start+window]
    score = sum(win) / len(win)
    numbers_to_plot.append(score)

plt.figure(figsize=(15,8))
plt.xlim(0,len(numbers_to_plot))
plt.plot(numbers_to_plot,linewidth=3,color="green")
plt.title('EXERCISE 2')
plt.ylabel('Unique amino acid residues')
plt.xlabel('Residue position')
plt.savefig("Chart_Exercise_2.png",transparent=True)
plt.show()




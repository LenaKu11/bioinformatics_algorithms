####################################################################################################
####################################################################################################
###########################################Assignment 2.############################################
####################################################################################################
####################################################################################################

################################################ A) ################################################

# Write a program that implements “patternCount(text, pattern)”, i.e. given a predetermined 
# pattern and a text, output how often the pattern occurs.

def patternCount(text, pattern):
    
    count = 0
    
    #Iterate over the string until the length of string - pattern
    for i in range (0, len(text) - len(pattern) + 1):
        
        #slice a sub-seq of the string
        sub_txt = text[i:i + len(pattern)]
        
        #compare the sub-seq with the pattern
        if sub_txt == pattern:
            count += 1
                
    return count
    
################################################ B) ################################################

# Count how many times the pattern 'ACGT' occurs

pwd = input('Input path to the file GCF_900086985.1_IMG-taxon_2623620618_annotated_assembly_genomic.fna in the form C:/user/.../ \n')

#Open the file reads.fna and put each line (=read) into a list
with open (f'{pwd}reads.fna', 'r') as read:
    reads = []
    for line in read:
        reads += [line.strip()]

def readsCount(text, pattern):
    #count how often the pattern occurs in a text which is stored in a list
    occ = 0
    
    for element in text:
        freq = patternCount(element, pattern)
        occ += freq
    return occ

#call the function readsCount for the given pattern and print out the result

pat = 'ACGT'
occ = readsCount(reads, pat)
all = len(reads)
pc = occ / all * 100

print(f'The pattern {pat} occures {occ} times in all {all} reads. This corresponds to {pc}% of the reads. ')

################################################ C) ################################################

# count in how many reads the pattern 'ACGT' occurs at least once

def countMin(sequence, pattern):
    occ = 0
    for element in sequence:
        freq = patternCount(element, pattern)
        if freq > 0:
            occ += 1
    return occ

pat = 'ACGT'   
occu = countMin(reads, pat)


print(f'The pattern {pat} occures in {occu} reads at least once.')


################################################ D) ################################################

#count the occurrence of all possible four base patterns

from itertools import product

#all possible combinations of('ACGT')
def possibleCombis(pattern):
    d = len(pattern)
    combinations = [''.join(i) for i in product(pattern, repeat = d)]
    return combinations

def CountOccurrence(text, pattern):
    counts = d = {}
    combis = possibleCombis(pattern)
    for element in combis:
        count = 0
        for line in text:
            freq = patternCount(line, element)
            count += int(freq)
        d = {element: count}
        counts.update(d)
    sort_comb = sorted(counts.items(), key=lambda x: x[1], reverse = True)

    return sort_comb
    
# call the functions for counting the occurences of all possible combinations of the four bases ACGT and print the result:
   
pat = 'ACGT'

mocc = CountOccurrence(reads, pat)[0:5]

print('The five most occuring patterns are: ')
for i in range (0, 5):
    print(*mocc[i], ' times')

################################################ E) ################################################

#output where the pattern occurs

from collections import Counter

def patternCount1(text, pattern):
    
    count = 0
    Index = []
    
    for element in text:
        
        ind = []
        for i in range (0, len(element) - len(pattern) + 1):
        
            sub_txt = element[i:i + len(pattern)]

            if sub_txt == pattern:
                count += 1
                ind.append(i)
                
        if len(ind) != 0:    
            Index.extend(ind)
    
    return count, Index


def indexGetter(text, pattern):
    
    patter = patternCount1(text, pattern)
    a = dict(Counter(patter[1]))
    index = sorted(a.items(), key=lambda x: x[1], reverse = True)
    
    return index

print('Frequency of the indices of the most abundant patterns in the reads: ')
for i in range (0, 5):
    b = mocc[i]
    c = indexGetter(reads, b[0])[0:5]
    print(f'{b[0]}: {c}')
####################################################################################################
####################################################################################################

################################################ A) ################################################

def HammingDistance(pattern1, pattern2):
    count = 0
    if len(pattern1) == len(pattern2):
        for char in range(0, len(pattern1)):
            if pattern1[char] != pattern2[char]:
                count += 1
            else:
                continue
        return count
            
            
def approximatePatternCount(text, pattern, d):
    #a function that counts how often a pattern occurs in a text with max d mismatches
    count = 0
    #ind = []
    for i in range(0, len(text) - len(pattern) + 1):
        seq = text[i:i + len(pattern)]
        if HammingDistance(seq, pattern) <= d:
            count += 1
            #ind.append(i - 1)
    return count#, ind



################################################ B) ################################################

pwd = input('Input path to the file GCF_900086985.1_IMG-taxon_2623620618_annotated_assembly_genomic.fna in the form C:/user/.../ \n')

#Open the file reads.fna and put each line (=read) into a list
with open (f'{pwd}reads.fna', 'r') as read:
    reads = []
    for line in read:
        reads += [line.strip()]

def appCountMin(sequence, pattern, d):
    occ = 0
    for element in sequence:
        freq = approximatePatternCount(element, pattern, d)
        if freq > 0:
            occ += 1
    return occ


pat = 'ACGT'
pat1= 'TGCA'
k = 1
n = 0

appCount1 = appCountMin(reads, pat, k) - appCountMin(reads, pat, n)
appCount2 = appCountMin(reads, pat1, k) - appCountMin(reads, pat1, n)
print(f'The pattern {pat} occures with {k} mismatch in {appCount1} reads at least once.')
print(f'The pattern {pat1} occures with {k} mismatch in {appCount2} reads at least once.')

################################################ C) ################################################

def errorRate(seq, pattern, n, k):
    
    appCountn = appCountMin(seq, pattern, n)    
    appCountk = appCountMin(seq, pattern, k) - appCountn
    eRate = appCountk / (appCountn + appCountk) * 100
    eRa = round(eRate, 2)
    
    return eRa

pat = 'ACGT'
err = errorRate(reads, pat, n, k)

print(f'The error rate for the pattern {pat} with one mismatch is {err} %')
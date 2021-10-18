####################################################################################################
####################################################################################################

# Write a program that calculates the GC-skew of a genome. Plot the GC-skew.                         #

from Bio import SeqIO
import matplotlib.pyplot as plt


pwd = input('Input path to the file GCF_900086985.1_IMG-taxon_2623620618_annotated_assembly_genomic.fna in the form C:/user/.../ \n')


for record in SeqIO.parse(f'{pwd}GCF_900086985.1_IMG-taxon_2623620618_annotated_assembly_genomic.fna', 'fasta'):
    halo = str(record.seq.upper())


def base_skew (sequence, B1, B2, step):
    
    count = 0
    skew_count = 0
    skew = []
    x = []
    
    for i in range (0, len(sequence), step):
        
        seq = halo[i:i + step]
        
        for char in seq:
            if char == B1:
                skew_count += 1
            
            elif char == B2:
                skew_count += - 1
            
            else:
                skew_count = skew_count
                
            count += 1
            
            x.append(count)
            skew.append(skew_count)
            
    y = f'{B1}{B2}_skew'
    plt.title(f'{B1}{B2} skew')   
    plt.xlabel('position')
    plt.ylabel('skew')
    plt.plot(x, skew)
    
    #show the coordinates of the extrema of the plot

    max_y = max(skew)
    min_y = min(skew)
    
    max_x = x[skew.index(max_y)]  # Find the x value corresponding to the maximum y value
    min_x = x[skew.index(min_y)]
    print ('max', max_y, 'max_x', max_x)
    print ('min', min_y, 'min_x', min_x)
    
    plt.show()
    
    plt.savefig(f'{pwd}{y}.png', format='png')


Ba1 = 'G'
Ba2 = 'C'
Ba3 = 'T'
wind = 1000

base_skew(halo, Ba1, Ba2, wind)

base_skew(halo, Ba1, Ba3, wind)

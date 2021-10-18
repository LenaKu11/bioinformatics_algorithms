##############################################################################
##############################################################################

# Solution to the Sectral Convolution Problem


def aaMass(aminoacid):
    
    aminoAcids = ['G', 'A', 'S', 'P', 'V', 'T', 'C', 'I', 'L', 'N', 'D', 'K', 'Q', 'E', 'M', 'H', 'F', 'R', 'Y', 'W']
    Masses = [57, 71, 87, 97, 99, 101, 103, 113, 113, 114,115, 128, 128, 129, 131, 137, 147, 156, 163, 186]
    aaMasses = dict(zip(aminoAcids, Masses))
    
    return aaMasses[aminoacid]

def singleMasses(peptide):
    
    spectrum = []
    for element in peptide:
        spectrum.append(aaMass(element))
        
    return spectrum
    

def PrefixMass(peptide):
    
    mass = [0]
    counter = 0
    masses = singleMasses(peptide)
    
    for i in range(len(masses)):
        counter += masses[i]
        mass.append(counter)
  
    return sorted(mass)


def cycloSpec(peptide):
    
    cycloSpectrum = [0]
    masses = PrefixMass(peptide) 
    peptideMass = masses[len(peptide)]
    
    for i in range(len(peptide)):

        for j in range(i + 1, len(peptide) + 1):
            cycloSpectrum.append(masses[j] - masses[i])
            if i > 0 and j < len(peptide):
                cycloSpectrum.append(peptideMass - (masses[j] - masses[i]))
    
    return sorted(cycloSpectrum)


##############################################################################

def convolutions(peptide):
    
    convolutions = []
    for i in range(len(peptide) - 1):
        temp = []
        for j in range(i + 1, len(peptide)):
            k = peptide[j] - peptide[i]
            if k > 0:
                temp.append(k)
        convolutions.append([peptide[i], sorted(temp)])
                
    return convolutions


def mostFreq(anestedList, spectrum):

    frequencies = []
    for i in range(len(anestedList)):
        for element in anestedList[i][1]:
            frequencies.append(element)
        frequencies.sort(reverse=True) 
     
    mostFreqEl = []
        
    for element in frequencies:
        if 57 <= element <= 200:
            mostFreqEl.append(element)
    mostFreqEl = sorted(set(mostFreqEl))
    
    frequencyCount = {}
    for i in range (len(mostFreqEl)):
        frequencyCount[mostFreqEl[i]] = 0
    
    for i in range(len(mostFreqEl)):
        count = 0
        for element in frequencies:
            if mostFreqEl[i] == element:
                count += 1
        frequencyCount[mostFreqEl[i]] = count
        
    return frequencyCount


def Mass():
    
    aminoAcids = ['G', 'A', 'S', 'P', 'V', 'T', 'C', 'I/L', 'N', 'D', 'K/Q', 'E', 'M', 'H', 'F', 'R', 'Y', 'W']
    Masses = [57, 71, 87, 97, 99, 101, 103, 113, 114,115, 128, 129, 131, 137, 147, 156, 163, 186]
    aaMasses = dict(zip(aminoAcids, Masses))
    
    return aaMasses



def translate(multiplicities):
    
    alist = []
    for element in multiplicities.items():
        alist.extend([[element[0], element[1]]])
        
    masses = Mass()  
    
    for element in masses.items():
        for i in range(len(alist)):
            if alist[i][0] == element[1]:
                alist[i] = ([element[0], alist[i][1]])

    return dict(alist)

        
        

def convolutionCycloPeptideSequencing(peptide):
    
    CycloSpectrum = peptide
    
    Convolutions = convolutions(CycloSpectrum)
    
    mostFrequentElements = mostFreq(Convolutions, CycloSpectrum)
    
    td = translate(mostFrequentElements)
    trl = sorted(td.items(), key=lambda item: item[1], reverse=True)
    
    #print('Most Frequent identified amino-acids are: ')
    #for i in trl:
    #    print(f'{i[0]} ({i[1]})')
    
    #mostFrequ = (sorted(mostFrequentElements.items(), key=lambda item: item[1], reverse=True))
    
    #print('Most frequent masses: ')
    #for j in mostFrequ:
    #    print(f'{j[0]} ({j[1]})')
        
    
    return dict(trl)

##############################################################################
##############################################################################
##############################################################################

theoretical_spectrum = [0, 57, 87, 97, 113, 144, 170, 184, 210, 241, 257, 267, 297, 354]
experimental_spectrum = [0, 57, 97, 99, 113, 144, 170, 184, 210, 241, 245, 257, 267, 297, 354]

print('multiplicities theoretical spectrum: ', convolutionCycloPeptideSequencing(theoretical_spectrum))
print('multiplicities experimental spectrum: ', convolutionCycloPeptideSequencing(experimental_spectrum))

##############################################################################
##############################################################################
##############################################################################
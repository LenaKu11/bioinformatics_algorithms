##############################################################################
##############################################################################

# Solution to LeaderBoardCycloPeptideSequencing from Bioinformatics Algorithms Textbook

def aaMass(aminoacid):
    
    aminoAcids = ['G', 'A', 'S', 'P', 'V', 'T', 'C', 'L', 'N', 'D', 'K', 'E', 'M', 'H', 'F', 'R', 'Y', 'W']
    Masses = [57, 71, 87, 97, 99, 101, 103, 113, 114, 115, 128, 129, 131, 137, 147, 156, 163, 186]
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

def linearSpec(peptide):
    
    linearSpectrum = [0]
    masses = PrefixMass(peptide) 
    
    for i in range(len(peptide)):
        for j in range(i + 1, len(peptide) + 1):
            linearSpectrum.append(masses[j] - masses[i])
            
    return sorted(linearSpectrum)

##############################################################################
##############################################################################
##############################################################################

#Expand: A collection containing all possible extensions of peptides in Peptides by a single amino acid mass
#set Peptides of length k, Expand construct set of length k + 1

def Expand(Board):
    
    #Here L means I/L, K means K/Q
    #aminoAcids = ['G', 'A', 'S', 'L', 'P', 'K']
    #aminoAcids = ['G', 'A', 'S', 'P', 'L', 'K', 'V', 'T']
    aminoAcids = ['G', 'A', 'S', 'P', 'V', 'T', 'C', 'L', 'N', 'D', 'K', 'E', 'M', 'H', 'F', 'R', 'Y', 'W']

    expanded = set()
    
    if Board == set('empty_peptide'):
        expanded = set(aminoAcids)
        
    else:
        for element in Board:
            for aminoacid in aminoAcids:
                expansion = element + str(aminoacid)
                expanded.add(expansion)

            
    return expanded


#Peptides total Mass/ Here: largest Mass in Spectrum:
def ParentMass(spectrum):
        return (sorted(spectrum, reverse=True)[0])

#Total Mass of the cyclic peptide
def Mass(peptide):
    return (sorted(cycloSpec(peptide), reverse=True)[0])
 

# Number of masses shared between Cyclospectrum (Peptide) and Spectrum (input)
# Take the multiplicities of shared masses in account 
def Score(peptide, spectrum):
    
    count = 0
    spec = list(spectrum)
    
    for value in peptide:
        for element in spec:
            if value == element:
                count += 1
                spec.remove(element)

    return count

def linearScore(peptide, Spec):
    peptidespec = linearSpec(peptide)
    return Score(peptidespec, Spec)
    

def Trim(leaderboard, Spectrum, N):
    
    leaderboard_scores = []
    LeaderBoard = set()
    
    for peptide in leaderboard:
        
        pep_spec = linearSpec(peptide)
        pep_score = Score(pep_spec, tuple(Spectrum))
        
        temp = []
        temp.append(peptide)
        temp.append(pep_score)
        leaderboard_scores.append(temp)  
    
    leaderboard_scores = sorted(leaderboard_scores, key = lambda x: x[1], reverse=True)
    
    if N >= len(leaderboard_scores):
        
        for i in range(len(leaderboard_scores)):
            LeaderBoard.add(leaderboard_scores[i][0])
            
    elif N < len(leaderboard_scores):
        
        for i in range(N):
            LeaderBoard.add(leaderboard_scores[i][0])  

            for j in range(N, len(leaderboard_scores)):
                if leaderboard_scores[j][1] == leaderboard_scores[N][1]:
                    LeaderBoard.add(leaderboard_scores[j][0])

    return LeaderBoard
        

def LeaderBoardCycloPeptideSequencing(Spectrum, N):
    
    leaderboard = set('empty_peptide')
    LeaderPeptide = ''
    mostFreq = []
    parentMass = ParentMass(Spectrum)
    
    while len(leaderboard) != 0:

        leaderboard = Expand(leaderboard)
        leaderboard_copy = list(leaderboard)

        for element in leaderboard:

            if Mass(element) == parentMass:
                
                element_spectrum = linearSpec(element)
                element_score = Score(element_spectrum, Spectrum)
                leaderpeptide_spectrum = linearSpec(LeaderPeptide)
                leaderpeptide_score = Score(leaderpeptide_spectrum, Spectrum)
                
                if element_score > leaderpeptide_score:
                    LeaderPeptide = element
                    temp = []
                    temp.append(LeaderPeptide)
                    temp.append(element_score)
                    mostFreq.append(temp)
                    
            elif Mass(element) > parentMass:
                leaderboard_copy.remove(element)
 
        leaderboard = Trim(leaderboard_copy, Spectrum, N)
        
    mostFrequentPeptides = sorted(mostFreq, key = lambda x: x[1], reverse=True)   
    print(mostFrequentPeptides)
    
    return LeaderPeptide

##############################################################################
##############################################################################
##############################################################################

t = [0, 57, 87, 97, 113, 144, 170, 184, 210, 241, 257, 267, 297, 354]
e = [0, 57, 97, 99, 113, 144, 170, 184, 210, 241, 245, 257, 267, 297, 354]

spectra = [t, e]
for element in spectra:
    print(element, ': ', LeaderBoardCycloPeptideSequencing(element, 1000))

##############################################################################
##############################################################################
##############################################################################
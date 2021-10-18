##############################################################################
##############################################################################

# Branch-and-bound algorithm for cyclopeptide sequencing and theoretical spectra
# write a code that outputs the linear and cyclic spectra of two peptides

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
        print()
        for j in range(i + 1, len(peptide) + 1):
            cycloSpectrum.append(masses[j] - masses[i])
            if i > 0 and j < len(peptide):
                cycloSpectrum.append(peptideMass - (masses[j] - masses[i]))
    
    return sorted(cycloSpectrum)
    
print('Spectrum of cyclic LFP: ', cycloSpec('LFP'))
print('Spectrum of cyclic VQY: ', cycloSpec('VQY'))

##############################################################################
##############################################################################
##############################################################################

def linearSpec(peptide):
    
    linearSpectrum = [0]
    masses = PrefixMass(peptide) 
    
    for i in range(len(peptide)):
        for j in range(i + 1, len(peptide) + 1):
            linearSpectrum.append(masses[j] - masses[i])
    
    return sorted(linearSpectrum)
    
print('Spectrum of linear LFP: ', linearSpec('LFP'))
print('Spectrum of linear VQY: ', linearSpec('VQY'))

##############################################################################
##############################################################################
##############################################################################
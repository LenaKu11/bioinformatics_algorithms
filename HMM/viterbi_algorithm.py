##############################################################################
##############################################################################

print('\n\nThis program generates the optimal state path and its probability \n(Viterbi Algorithm), according to the transition and emission probabilities given in the Algorithmic Bioinformatics lecture. ')
inp = input('Insert your sequence (H/T) here and press enter: ')
print('\n\n...calculating state path...')

def translated(obs, OSpace):
    #decode heads and tails to 0 and 1, since those can be assigned to matrix entries
    obsDecoded = []
    for element in obs:
        if element == OSpace[0]:
            obsDecoded.append(0)
        else:
            obsDecoded.append(1)

    return obsDecoded

def dptable(VMat):
    # Print a table of steps from dictionary
    yield " " * 5 + "     ".join(("%3d" % i) for i in range(len(VMat)))
    for state in VMat[0]:
        yield "%.7s: " % state + " ".join("%.7s" % ("%lf" % VMat[state]["pr"]) for VMat in VMat)

def viterbi(OSpace, SSpace, observations, states, probs, Tm, Em):
    obs = translated(observations, OSpace)

    VMat = [{}]
    for state in states:
        VMat[0][state] = {"pr": probs[state] * Em[state][obs[0]], "pre": None}
    
    # Run Viterbi when t > 0
    for t in range(1, len(obs)):
        VMat.append({})
        for state in states:
            max_tr_prob = VMat[t - 1][states[0]]["pr"] * Tm[states[0]][state]
            prev_st_selected = states[0]
            for prev_st in states[1:]:
                tr_prob = VMat[t - 1][prev_st]["pr"] * Tm[prev_st][state]
                if tr_prob > max_tr_prob:
                    max_tr_prob = tr_prob
                    prev_st_selected = prev_st

            max_prob = max_tr_prob * Em[state][obs[t]]
            VMat[t][state] = {"pr": max_prob, "pre": prev_st_selected}
    print('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
    print('input sequence: \n\n', observations)
    #Tprint('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')  


    #for line in dptable(VMat):
    #    print(line)

    opt = []
    max_prob = 0.0
    best_st = None
    
    # Get most probable state and its backtrack
    for state, data in VMat[-1].items():
        if data["pr"] > max_prob:
            max_prob = data["pr"]
            best_st = state
    opt.append(best_st)
    previous = best_st

    # Follow the backtrack till the first observation
    for t in range(len(VMat) - 2, -1, -1):
        opt.insert(0, VMat[t + 1][previous]["pre"])
        previous = VMat[t + 1][previous]["pre"]
    
    #re-translate the sequence from 0/1 to F/B
    o = ''
    for element in observations:
        o += element

    tr = ''
    for element in opt:
        if element == 0:
            tr += SSpace[0]
        else:
            tr += SSpace[1]
            
    #print('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
    #print('input sequence: \n\n', observations)
    print('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')        
    print('optimal state path for the sequence: \n\n', tr, '\n\nwith probability ', max_prob)
    print('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    
    return
        
O = ['H', 'T']     # O: Observation space
S = ['F', 'B']     # S: State space

probs = [0.5, 0.5] # probs: initial probabilities for the states F/B

Tm = [[0.9, 0.1], 
      [0.1, 0.9]]  # Tm: Transmission matrix

Em = [[0.5, 0.5], 
      [0.75, 0.25]]# Em: Emission Matrix


state = [0, 1]

#emSeq = ['THHTTTTHHTHHHHTTHHHH', 'HTTHHHTHHHTHHTHHTHTH', 'TTHHTTHHTHTHHHHHHHTT']

#for element in emSeq:
#    print(viterbi(O, S, element, state, probs, Tm, Em))

viterbi(O, S, inp, state, probs, Tm, Em)
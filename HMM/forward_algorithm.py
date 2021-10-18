##############################################################################
##############################################################################

print('This program computes the likelihood of the HMM given in the Bioinformatics Lecture to emit a given sequence.')
observations = input('Input a sequence (H/T): ')


def forward(observations, states, init, Tm, Em, end_st):
    fwd = []
    for i, observation_i in enumerate(observations):
        currF = {}
        for st in states:
            if i == 0:
                preSum = init[st]
            else:
                preSum = sum(preF[k] * Tm[k][st] for k in states)

            currF[st] = Em[st][observation_i] * preSum

        fwd.append(currF)
        preF = currF

    p_fwd = sum(currF[k] * Tm[k][end_st] for k in states)
    
    print('\nInput sequence: \n', observations)
    print('\nLikelihood of emitting a symbol for a state at each position of the sequence: ')
    
    print('F\t\t\t B')
    for element in fwd:
        print(round(element['F'], 4), '\t\t', round(element['B'], 4))
    
    print('\nProbability of the given HMM to emit this sequence:\n', p_fwd)
    
    return fwd, p_fwd

states = ('F', 'B')
end_state = 'E'
 
start_probability = {'F': 0.5, 'B': 0.5}
 
transition_probability = {'F' : {'F': 0.9, 'B': 0.1, 'E': 0.01}, 'B' : {'F': 0.1, 'B': 0.9, 'E': 0.01},}
 
emission_probability = {'F' : {'H': 0.5, 'T': 0.5},'B' : {'H': 0.75, 'T': 0.25},}
#emSeq = ['THHTTTTHHTHHHHTTHHHH', 'HTTHHHTHHHTHHTHHTHTH', 'TTHHTTHHTHTHHHHHHHTT']
#for observations in emSeq:

forward(observations, states, start_probability, transition_probability, emission_probability, end_state)
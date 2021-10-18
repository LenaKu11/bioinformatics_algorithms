##############################################################################
##############################################################################

print('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
print('\nThis program creates emission sequences (Heads/Tails) and state paths \n(Fair/Biased) according to the transition and emission probabilities\ngiven in the Algorithmic Bioinformatics lecture. ')
print('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
inp = input('\nEnter desired sequence length: ')
numb = input('\nHow many sequences do you need? ')
print('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')

import numpy as np

def emission(state):
    states = ['H', 'T']
    if state == 'F':
        probs = [0.5, 0.5]
    if state == 'B':
        probs = [0.75, 0.25]
    return np.random.choice(states, p=probs)

def transition(state):
    states = ['F', 'B']
    if state == 'F':
        probs = [0.9, 0.1]
    if state == 'B':
        probs = [0.1, 0.9]
    return np.random.choice(states, p=probs)

def sequenceGen(k):
    state = 'F'
    statePath = state
    emissionSeq = ''
    for i in range(k):
        emissionSeq += emission(state)
        state = transition(state)
        statePath += state
        
    print('\nEmitted Sequence:\n', emissionSeq)
    print('\nState Path:\n', statePath[:-1])

    
    return emissionSeq, statePath[:-1]


for i in range(int(numb)): 
    print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    print(f'\nSequence number {i + 1}:')
    sequenceGen(int(inp))
    print()
    
print('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
##############################################################################


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
    
    return fwd, p_fwd


def backward(observations, states, init, Tm, Em, end_st):
    
    back = []
    
    for i, observation_i_plus in enumerate(reversed(observations[1:] + '0')):
        currB = {}
        for st in states:
            if i == 0:
                currB[st] = Tm[st][end_st]
            else:
                currB[st] = sum(Tm[st][l] * Em[l][observation_i_plus] * preB[l] for l in states)

        back.insert(0,currB)
        preB = currB

    p_bkw = sum(init[l] * Em[l][observations[0]] * currB[l] for l in states)

    return back


def merge(observations, states, start_probability, transition_probability, emission_probability, end_state):
    
    fwd_start = forward(observations, states, start_probability, transition_probability, emission_probability, end_state)
    fwd = fwd_start[0]
    p_fwd = fwd_start[1]
    bkw = backward(observations, states, start_probability, transition_probability, emission_probability, end_state)
    
    posterior = []
    for i in range(len(observations)):
        posterior.append({state: fwd[i][state] * bkw[i][state] / p_fwd for state in states})
    
    print('\nInput sequence: \n', observations)
    print('\nResponsibility Matrix:')
    print()
    print('  |\tF \t\t B')
    print('------------------')
    for i, element in enumerate(posterior):
        print(observations[i], '|', round(element['F'], 3), '\t', round(element['B'], 3))
           
    return #posterior

inp = input('Input a emission sequence (H/T): \n')

states = ('F', 'B')
end_state = 'E'
 

 
start_probability = {'F': 0.5, 'B': 0.5}
 
transition_probability = {'F' : {'F': 0.9, 'B': 0.1, 'E': 0.01}, 'B' : {'F': 0.1, 'B': 0.9, 'E': 0.01},}
 
emission_probability = {'F' : {'H': 0.5, 'T': 0.5},'B' : {'H': 0.75, 'T': 0.25},}


merge(inp, states, start_probability, transition_probability, emission_probability, end_state)
    
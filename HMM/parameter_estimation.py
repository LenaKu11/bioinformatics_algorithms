##############################################################################
##############################################################################

print('\n\nParameter Estimation for HMMs.\n\n')
sequence = input('Input Sequence of symbols (H or T): ')
statepath = input('Input State Path of symbols (F or B): ')

def EmCalc(seq):
    length = len(seq)
    
    heads = 0
    tails = 0
    for i in range(length):
        if seq[i] == 'H':
            heads += 1
        elif seq[i]:
            tails += 1
    
    if length == 0:
        h, t = 0.0, 0.0
        
    else:
        h, t = round(heads/length, 2), round(tails/length, 2)
    
    return [h, t]

def TmCalc(seq, symbol):
    
    if symbol not in seq:
        return [0.0, 0.0]
    
    length = len(seq)
    cont = 0
    discont = 0

    for i in range(0, length - 1):
        if seq[i] == symbol:
            if seq[i + 1] == symbol:
                
                cont += 1
            else:
                discont += 1

    l = (cont + discont)
    if symbol == 'F':
        c, d = round(cont/l, 2), round(discont/l, 2)
    else:
        d, c = round(cont/l, 2), round(discont/l, 2)

    
    return [c, d]

def dev(param, state, ind, Mat):
    result = 0
    if Mat == 'Tm':
        if state == 'F':
            if ind == 0:
                result = (param - 0.9)
            else:
                result = (param - 0.1)
        else:
            if ind == 1:
                result = (param - 0.9)
            else:
                result = (param - 0.1)
                
    if Mat == 'Em':
        if state == 'F':
            result = (param - 0.5)
        else:
            if ind == 0:
                result = (param - 0.75)
            else:
                result = (param - 0.25)
    
    res = round(result, 2)
    
    if res >= 0:
        res = f'+{res}'
                
    return res

def paramEstimation(sequence, statePath):
    
    if len(sequence) != len(statePath):
        print('Sequence and State Path must have equal length.')
        return
    
    print('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    print('Input Sequence: \n', sequence)
    print('Input State Path:\n', statePath)
    length = len(sequence)
    fair = []
    biased = []
    for i in range(length):
        if statePath[i] == 'F':
            fair.append(sequence[i])
        else:
            biased.append(sequence[i])
    
    #Calculate Emission probability:
    f = EmCalc(fair)
    b = EmCalc(biased)
    
    Em = [f, b]

    #calculate Transition probabilities
    
    f = TmCalc(statePath, 'F')
    b = TmCalc(statePath, 'B')
    
    Tm = [f, b]

    print('\n~~~~~~~~~~~~~~~~~~~~~Estimated parameters:~~~~~~~~~~~~~~~~~~~~~ \n')
    print('(and the deviation from the parameters given in the lecture, in brackets)\n\n')
    print('Emission Matrix: ', '\n')
    
    print('\t H \t\t\t\t T')
    
    states = ['F', 'B']
    for i in range(len(Em)):
        print(states[i], '\t', Em[i][0], '(', dev(Em[i][0], states[i], 0, 'Em'), ')\t', Em[i][1], '(', dev(Em[i][1], states[i], 1, 'Em'), ')')

    print('\n Transition Matrix: ', '\n')
    print('\t F \t\t\t\t B')
    for i in range(len(Tm)):
        print(states[i], '\t', Tm[i][0], '(', dev(Tm[i][0], states[i], 0, 'Tm'), ')\t', Tm[i][1], '(', dev(Tm[i][1], states[i], 1, 'Tm'), ')')
    print('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    
    return Em, Tm
                
sequence = 'HHHHTHHTTHHTTTTTTTTHHHTTTHHTHHHHHTHTTHHTHHTHHTTHTHTTTTTHTHHHHHTHHTTHTTTHTHTHHHHHHHHHHHTHTHTTHHHTHHTT'

statepath = 'FFFBBBFFFFFFFFFFFFFFFFFFFFFFFFBBBFFBBBBBBBBBBBBFFFFFFFFFFFFFFFFFFFFFFFFFFFFBBBBBBBBBBBBBBBBFFFFFFFFF'
    
paramEstimation(sequence, statepath)

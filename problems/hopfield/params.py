import scipy

inputState = [1,1,-1,-1,1]
# inputState = [1,-1,1,-1]
# inputState = [1,-1,1]
numQubits = len(inputState)
includeInput = True
numMemories = 4
simCase = str(numMemories)
numEigvals = 5
#annealTime = 3.0
annealTime = scipy.arange(0.1,15,0.5)

learningRules = ['hebb', 'stork', 'proj']
rule = 2

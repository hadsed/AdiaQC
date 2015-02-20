'''
Generate some Hopfield problem using Hebb rule.
'''

#TODO: embed this in the config generator
import numpy as np

# hopfield params
neurons = 8
# first half 0's, second half 1's
#inputstate = [1]*(neurons/2) + [-1]*(neurons/2)
# one memory that's all ones
#could go up to 4 mems with n=32?
memories = [ [1]*neurons ]
#memories = [ [1]*neurons, [1,-1]*(neurons/2), [1,1,-1,-1]*(neurons/4), [1]*(neurons/2) + [-1]*(neurons/2-1)+[1] ]
print memories
# hebb rule
memMat = np.matrix(memories).T
ising_off1 = (np.triu(memMat*memMat.T) - 
             len(memories)*np.eye(neurons))/float(neurons)
ising_off2 = np.triu(memMat*memMat.T, 1)/float(neurons) 
#ising_diag = np.array(inputstate)
print ising_off1
print
print ising_off2
print
print ising_off2[2,3]
#print ising_diag
# save to file
#np.save('hopfield'+str(neurons)+'_mat.npy', ising_off)
#np.save('hopfield'+str(neurons)+'_inp.npy', ising_diag)
# test to see that it works right
#io = np.load('hopfield'+str(neurons)+'_mat.npy')
#id = np.load('hopfield'+str(neurons)+'_inp.npy')
#print("Saved file same as generated?")
#print(np.all(io == ising_off) and np.all(id == ising_diag))

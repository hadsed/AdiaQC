'''
Walk through all instance directories and run a script in there.
'''

import os

cmd = 'python ../../../../tools/hopfield_analysis.py -p probsT15.0.dat -b 0 -k 2'

for root, dirs, files in os.walk('.'):
    if root == '.':
        continue
    if dirs == []:
        os.chdir(root)
        os.system(cmd)
        os.chdir('../../')

AdiaQC
======

Adiabatic quantum computation solver used to solve various problems (machine learning problems included). This is written in Python 3.2 and tested on Linux only, but may work for you (do let me know if it does). Requires SciPy/NumPy.

Running a problem looks like:

python3 run.py -p [problem]

where [problem] can be hopfield.py, which is included. Outputs final probabilities with qubit labelings included. You should only edit the problem files in problems/, nothing else should be changed unless you know what you're doing.

Problems are included in the problems/ directory. Currently working on Hopfield neural network, correlation clustering to come soon. Add your own problem by editing problems/problems.py and adding an initialization function that outputs in the same format as the others. Also, add a new file that includes similar parameters as in problems/hopfield.py. 'nQubits', 'Q', and 'a' must all be defined in terms of the problem--how that happens does not matter.

One can also add a new solver by creating a new function in solver.py and calling it from run.py.

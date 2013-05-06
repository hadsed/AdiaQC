AdiaQC
======

Adiabatic quantum computation solver used to simulate various problems (some machine learning problems included). This is written in Python 3.2 and tested on Linux only, but may work for you (do let me know if it does). I try to keep it compatible with Python 2.6+, but not actively. Requires SciPy (NumPy, PyLab).

Running a problem looks like:

python3 run.py -p [problem]

where [problem] can be hopfield.py, which is included (the *.py is redundant, the extension is removed anyway). Outputs final probabilities with qubit labelings included. You should only edit the problem files in problems/, nothing else should be changed unless you're trying to modify the code.

Please run the test problems. They are simple problems that you can verify by hand or by a brute-force QUBO solver. Some tell what the answers should be in the comments.

Problems are included in the problems/ directory. Currently working on Hopfield neural network, correlation clustering to come soon. Add your own problem by creating a new problem in the problems/ directory. 'nQubits', 'Q', 'a', 'dt', and 'T' must all be defined in terms of the problem. All of those parameters are floats except for 'Q', which is the upper-triangular QUBO matrix. Output parameters like 'eigspecflag' and 'outputdir' are required too.

You don't have to worry about namespace cluttering, the only variables taken from your problem script will be the ones that are required. Feel free to do anything you want (import packages, read from files, etc.) to help you get those parameters inside your problem file.

You can also add a new solver by creating a new function in solver.py and calling it from run.py.

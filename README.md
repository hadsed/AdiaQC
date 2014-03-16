AdiaQC
======

Adiabatic quantum computation solver used to simulate various problems (some machine learning problems included). This is written in Python 2.7.4 and tested on Linux only, but may work for you (do let me know if it does). Requires SciPy (NumPy, PyLab). Probably works with Python3 with some little tweaks (if you make them, please share).

Running a problem looks like:

python2 run.py -p [problem]

where [problem] can be hopfield.py, which is included (the *.py is redundant, the extension is removed anyway and you can specify without it).

The problem files define a function called 'parameters'. It takes in command-line arguments as a dictionary. This is useful for simulations with large numbers of instances if you need to do something special for each instance, or you want to specify a simulation type. The cmd args are:

* --problem: specify the problem file path relative to AdiaQC/problems
* --relpath: this should be specified so that output files are put in the right place if you're not running run.py from AdiaQC/
* --instance: you can specify an instance number, if it matters, e.g. with the hopfield network, if we want to try all possible combinations of P-memory networks, we can generate a list (a script for this is included in tools/) and then assign each instance an item in this list so that you can specify all the instances you want without having to make a thousand problem files separately for each one
* --simtype: this is another similar parameter to --instance. it's useful, for example, if you have a small modification to how you build your problem, like in hopfield networks with different learning rules. then in the problem file you can write if-else statements to do the appropriate thing. again, this means that you don't have to create extra problem files and you can keep it all in one place.

You should only have to edit the problem files in problems/, nothing else should be changed unless you're trying to modify for your own purposes.

Please run the test problems. They are simple problems that you can verify by hand or by a brute-force QUBO solver (though the answers should be at the top if they are simple). If you have cool problems, feel free to contribute them.

Problems are included in the problems/ directory. Currently working on Hopfield neural network, correlation clustering to come soon. Add your own problem by creating a new problem in the problems/ directory. At the end of any problem file is a list of what is required to be returned by that function. All of them need to be specified. Feel free to do anything you want inside that function (import packages, read from files, etc.) to help you get those parameters.

If you are looking to do large numbers of instances, do look at what I did with the Hopfield network stuff. The relevant files are problems/hopfield_batch.py and problems/hopfield_wrapper.py. The latter is the wrapper script which runs a specified instance of hopfield_batch as its own process.

You can also add a new solver by creating a new function in solver.py and calling it from run.py. If you do that, please submit a pull-request and contribute it!

If you are looking to use this code or modify it for your purposes, please don't hesitate to contact me through email. My address is the first 3 letters of my first name followed by the first 3 letters of my last name at google mail.

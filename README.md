# IE
C code to solve integral equations for primordial density perturbations


These are the codes I used to solve the evolution of primordial perturbations as described in arXiv:2105:?????.

There are two codes:  iterative.c and boltzmann_it.c.  They both evolve primordial perturbations until after recombination but before reionization.  I am making these codes avaiable so anyone who is interested in details of the calculation described in the paper can find them.  The code is *not* certified for cosmological analysis.  The code uses a (poor implentation of) the UFA approximation; and uses an analytic approximation for the expansion rate at early times.  It uses an ionization fraction calculated from from HyRec-.

In addition to iterative.c and boltzmann_it.c, there are a number of include files.  

The output from HyRec-2 is in output_xe.dat.  The input from a prior run of either code is in iterative.in.

I compiled both codes on my Mac laptop using

gcc -g iterative.c -lm -o iterative
gcc -g boltzmann_it.c -lm -o boltzmann_it

although I think that my gcc is actually clang


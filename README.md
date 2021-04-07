# IE
C code to solve integral equations for primordial density perturbations


This is the code I used to solve the evolution of primordial perturbations as described in arXiv:21??:?????.

The code evolves primordial perturbations until after recombination but before reionization.  The code is provided as is so interested parties can find details of the calculation described in the paper.  The code is *not* certified for cosmological analysis.  The code sets the baryon sound speed to zero; uses a (poor implentation of) the UFA approximation; and uses an analytic approximation for the expansion rate at early times.  It uses an ionization fraction calculated from elsewhere---I used HyRec-2.  The numerics also need to be developed and checked further.

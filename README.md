# phasta_SST
SST k-w model implementation for PHASTA - CFD course final project.

All the e3 routines along with getdiff.f are to be added to phSolver/incompressible.

input_fform.cc, turbKW.f, gradVgen.f and elm3komega are to be added to phSolver/common.

Careful while using a solver.inp file. (input_fform.cc has been modified to use user-defined values 
of characteristic velocity and freestream lengths. These two features are, without implementation
of boundary conditions, currently redundant.)

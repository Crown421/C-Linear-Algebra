# Scientific Computing with C++ #

This repository contains the code for my Scientific Computing project at the university of Oxford.
The task was comprised of two tasks, firstly implementing a matrix class and relevant
methods, and secondly an extension of our own choice. I chose pursue various optimization
approaches and compared them to the performance of Matlab code.

## The matrix class ##

For matrix storage I chose dynamic allocation of arrays of arrays containing the rows, so
row-based storage scheme as sketched in the following subfigure **A**. Taking inspiration from the Matlab syntax,
I further implemented slicing operators using round brackets. The facilitate both the efficient use of storage and
the change of submatrices of a matrix, shallow copies where used. These are matrix instances consisting just of pointers and arrays of pointers, referencing the original content, as also visible in the figure. This potentially dangerous, so various safeguards where implemented to avoid corruption of data and conflict with return value optimization.

<!-- ![Storage method](https://crown421.github.io/rep_hosting/C-Linear_Algebra/class.png) -->
<img src="https://crown421.github.io/rep_hosting/C-Linear_Algebra/class.png" width="750">

Further, the main algorithm in this work, GMRES, requires matrices of changing size for each step. As part of the optimization, an alternative constructor was written to allocate a large matrix.
As sketched in the previous subfigure **B**, only a part of the allocated storage may be accessible to arithmetic and element-altering methods.  

## Algorithm: LU decomposition

The first implemented solver was the [LU decomposition](https://en.wikipedia.org/wiki/LU_decomposition), which was implemented as a separate class. The constructor takes a matrix object, and returns an object containing the LU decomposition. The class methods return various related things, like the solution for a given matrix of right-hand side columns or the matrix's eigenvalue.

The LU class is explicitly not a friend, but works with the matrix class's getters/setters and its arithmetic methods. Those were implemented using safeguards like bounds checking, which are not necessary once the class is functional.
Hence, it is possible to disable those safeguards via a compiler flag, which is indicated by the suffix ``S`` in the following timing plots.
Additional optimization was done by additionally compiling with ``O3``.

The results are timed for a purely random matrix with positive entries between ``0`` and ``1``, the same matrices with ``5`` added to the diagonal, and lastly matrices with only 5 distinct eigenvalues. For comparison, the build-in Matlab solver ``mldivide`` is also included. For each size, 10 different matrices are solved 10 times and the results averaged.

<!-- ![LU performance](https://crown421.github.io/rep_hosting/C-Linear_Algebra/LUsolve.png) -->
<img src="https://crown421.github.io/rep_hosting/C-Linear_Algebra/LUsolve.png" width="750">

## Algorithm: GMRES
The [GMRES algorithm](https://en.wikipedia.org/wiki/Generalized_minimal_residual_method) is a sophisticated algorithm, best suited for general, very large, sparse matrices with clustered eigenvalues. It is a Krylov-space methods, which means that for each step another basis vector for the solution space is computed.
Hence, two versions where implemented. In ``1`` the dynamically allocated space is resized and hence recopied in each step, and in ``2`` the maximal necessary space is reserved in the beginning.

Those algorithms are timed on the same matrices, using the same optimization as in the previous section. They are further compared to the Matlab build-in function ``gmres`` (BIM) and most importantly an implementation of the same algorithm using Matlab (MSL).
<!-- ![GMRES performance](https://crown421.github.io/rep_hosting/C-Linear_Algebra/gmres.png) -->
<img src="https://crown421.github.io/rep_hosting/C-Linear_Algebra/gmres.png" width="750">

# DQGMRES
**D**irect **Q**uasi-minimal **GMRES**  
<sup>_#numericallinearalgebra_ | _#krylovmethods_</sup>
___
### Fast Krylov subspace method for indefinite square systems

Use of incomplete orthogonalization[^1] leads to reduced memory demands compared to the 'full' GMRES algorithm, which stores the entire Krylov Matrix.
Adapted from GMRES implemetation in [IterativeSolvers.jl](https://github.com/JuliaMath/IterativeSolvers.jl) (this package is in fact a dependency).


[^1]: Saad, Y. and Wu, K., 1996. DQGMRES: A direct quasi‐minimal residual algorithm based on incomplete orthogonalization. Numerical linear algebra with applications, 3(4), pp.329-343.

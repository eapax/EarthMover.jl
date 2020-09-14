[![Build Status](https://travis-ci.com/eapax/EarthMover.jl.svg?branch=master)](https://travis-ci.com/eapax/EarthMover.jl)

# EarthMover.jl

This is a package for computing optimal transport distances, such as the Wasserstein distances, between probability measures.

These distances are motivated by thinking of probability measures as mass distributions. The distance between two probability measures p and q is defined as the lowest cost at which one can transport mass from p to q---having introduced a *cost function* c(x,y) which defines the cost to move unit mass from x to y. The case c(x,y)=|x-y|^p defines the p-Wasserstein distance.

These distances are the most natural way to compare between probability measures. Indeed, a famous theorem says that the Wasserstein distance *metrizes* the space of probability measures ([1, Theorem 7.3]). But computation of these distances leads to challenges.

The problem of Optimal Transport has two alternative formulations: that of Monge and that of Kantorovich. The Monge formulation treats discrete probability distributions as collections of N equal masses (think of N shipping crates on a warehouse floor), and defines a possible transport plan as a permutation of the N objects. The Kantorovich formulation is a relaxation of this which applies to continuous probability distributions (think piles of sand), and allows for transport plans which subdivide the mass arbitrarily. The two formulations turn out to give the same result (i.e. the same definition of "distance") whenever both are defined ([1, pg5]), and both formulations lead to difficult (though spiritually different!) optimization problems.

The Monge formulation leads one to a famous problem in combinatorics, the *assignment problem*, whilst the Kantorovich formulation leads one to a classic problem in linear optimization. The former can be solved (in polynomial time) by the famous Hungarian algorithm, whilst the latter can be solved (usually in polynomial time) by one of the many well-developed linear programming software packages. In practice, one will find that one method is more effective than another for the problem at hand. 

***This package provides the functionality for both methods, and more.***

Specifically, this package exports functions to compute optimal transport distances by implementing one of the following methods "under-the-hood":

1) Solution of the assignment problem by the Hungarian algorithm (via the native Julia dependancy Hungarian.jl)
2) Solution of the linear optimization problem (via the well-optimized Python optimization package scipy.optimize.linprog, called via the dependency PyCall.jl)
3) Approximate solution of the linear optimization problem by the Sinkhorn algorithm as described in [2] (implemented here in native Julia, without any dependencies).

Note that 1 and 2 give *exact* computation of optimal transport distances, whilst 3 gives an *approximate* computation (but is often much faster than 1 or 2).

## Usage

A simple usage example coming soon.


## Installation

EarthMover.jl is not registered in the Julia Registry, so do
```julia
julia>] add https://github.com/eapax/EarthMover.jl
```
where `]` opens the package manager.

## References

[1]: C. Villani -- Topics in optimal transport. 

[2]: M. Cuturi -- Sinkhorn distances: lightspeed computation of optimal transport distances. 

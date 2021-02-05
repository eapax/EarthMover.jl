[![Build Status](https://travis-ci.com/eapax/EarthMover.jl.svg?branch=master)](https://travis-ci.com/eapax/EarthMover.jl)

# EarthMover.jl

A package for computing optimal transport (OT) distances, such as the Wasserstein distances, between probability measures.

These distances are motivated by thinking of probability measures as mass distributions. The distance between two probability measures p and q is defined as the lowest cost at which one can transport mass from distribution p to distribution q---having introduced a *cost function* c(x,y) which defines the cost to move unit mass from x to y. The case c(x,y)=|x-y|^p defines the p-Wasserstein distance.

These distances are the most natural way to compare between probability measures. Indeed, a famous theorem says that the Wasserstein distance *metrizes* the space of probability measures ([1, Theorem 7.3]). But computation of these distances leads to challenges.

This package exports functions to compute OT distances via one of the following methods:

1) Solution of the assignment problem (cf. the Monge formulation of OT) by the Hungarian algorithm (via the native Julia Hungarian.jl)
2) Solution of the linear optimization problem (cf. the Kantorovich formulation of OT) by calling the Python optimization function scipy.optimize.linprog (via PyCall.jl)
3) Approximate solution of the linear optimization problem by the Sinkhorn algorithm as described in [2] (implemented here in native Julia, without any dependencies).

Note that 1 and 2 give *exact* computation of optimal transport distances, whilst 3 gives an *approximate* computation (often much faster than 1 or 2). Which method, 1, 2 or 3, is most appropriate will depend on the problem at hand. 

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

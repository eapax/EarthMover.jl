[![Build Status](https://travis-ci.com/eapax/EarthMover.jl.svg?branch=master)](https://travis-ci.com/eapax/EarthMover.jl)

# EarthMover.jl

A package for computing optimal transport distances, such as the Wasserstein distances (WDs), between probability measures.

These distances are motivated by thinking of probability measures as mass distributions. The distance between two probability measures p and q is defined as the lowest cost at which one can transport mass from distribution p to distribution q---having introduced a *cost function* c(x,y) which defines the cost to move unit mass from x to y. The case c(x,y)=|x-y|^p defines the p-Wasserstein distance.

These distances are the most natural way to compare between probability measures. Indeed, a famous theorem says that the Wasserstein distance *metrizes* the space of probability measures ([1, Theorem 7.3]). But computation of these distances leads to challenges.

This package exports functions to compute distances via one of the following approaches:

1) Solution of the assignment problem (cf. the Monge formulation [1]) by the Hungarian algorithm (via the native Julia Hungarian.jl)
2) Solution of the linear optimization problem (cf. the Kantorovich formulation [1]) by calling the Python optimization function scipy.optimize.linprog (via PyCall.jl)
3) Approximate solution of the linear optimization problem by the Sinkhorn algorithm as described in [2] (implemented here in native Julia, without any dependencies).

Note that 1 and 2 give exact computation of WDs, whilst 3 computes Sinkhorn divergences (SDs) which approximate WDs (often faster than 1 or 2). Which method, 1, 2 or 3, is most appropriate will depend on the problem at hand. 

Currently this package is only implemented with the p=1 case c(x,y)=|x-y| though it would be simple to add other costs.

Also note that for one dimensional data sets there exists an explicit formula for the WD, 
which can be computed extremely fast, but this is not yet implemented here.
Thus for 1D data it is much better to 
use an alternative package such as `scipy.stats.wasserstein_distance` for Python.

## Usage

Approach 1 can be called through `solve_monge` on a pair of raw data sets of the same length, computing the WD between the corresponding empirical distributions. For example, using some 2D data

```julia
julia> X = [(2.0, 3.0), (4.0, 5.0), (6.0, 1.0)]
julia> Y = [(1.0, 3.0), (4.0, 6.0), (7.0, 1.0)]
julia> solve_monge(X, Y)
```
will return `1.0` which is easily checked by hand. The data can alternatively be binned, before passing the histograms to the Kantorovich formulation (approach 2)

```julia
julia> using StatsBase
julia> edges = (0:10, 0:10)
julia> Xdata = ([2.0, 4.0, 6.0], [3.0, 5.0, 1.0])
julia> Ydata = ([1.0, 4.0, 7.0], [3.0, 6.0, 1.0])
julia> H1 = fit(Histogram, Xdata, edges)
julia> H2 = fit(Histogram, Ydata, edges)
julia> solve_kantorovich(H1, H2)
```

This approach is useful for large data arrays as the complexity depends only on the binwidth. It is also applicable when X and Y are data sets of different lengths, where the Monge formulation does not apply.

The Sinkhorn divergence (SD) (approach 3) can be called via

```julia
julia> epsilon = 0.05
julia> solve_sinkhorn(X, Y, epsilon)
julia> solve_sinkhorn(H1, H2, epsilon)
```

Note that this works either with raw data (X,Y) or binned data (H1,H2). 
For sufficiently small epsilon the SD will approximate the WD, 
and may provide significant speed-ups for epsilon not too small [2].

Currently `solve_kantorovich` and `solve_sinkhorn` only accept histograms of dimension 1, 2, or 3. 
This is somewhat intentional since for high-dimensional data-sets binning becomes inefficient, 
so working with raw unbinned data through `solve_monge` or `solve_sinkhorn` is preferred.
Note, however, that the WD suffers from a "curse of dimensionality" in this context, 
while the SD may not [3].


## Installation

EarthMover.jl is not registered in the Julia Registry, so do
```julia
julia>] add https://github.com/eapax/EarthMover.jl
```
where `]` opens the package manager.

## References

[1]: C. Villani -- Topics in optimal transport (2003).

[2]: M. Cuturi -- Sinkhorn distances: lightspeed computation of optimal transport distances (2013). 

[3]: A. Genevay, L. Chizat, F. Bach, M. Cuturi, G. Peyré -- Sample Complexity of Sinkhorn Divergences (2019).

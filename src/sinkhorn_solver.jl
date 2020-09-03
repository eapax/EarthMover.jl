#This code implements the "Sinkhorn algorithm" as it is described in the wonderful paper of Cuturi: "Sinkhorn Distances: Lightspeed Computation of Optimal Transport"

using StatsBase

include("histo_unwrap_kantorovich.jl")

function solve_sinkhorn(H1::Histogram, H2::Histogram, eps::Real=0.125)
    @assert eps>0
    @assert length(H1.edges)==length(H2.edges)<=3
    m = length(H1.edges)
    if m==1
        (P,Q,c) = make_probcost_vectors1d(H1,H2)
    end
    if m==2
        (P,Q,c) = make_probcost_vectors2d(H1,H2)
    end
    if m==3
        (P,Q,c) = make_probcost_vectors3d(H1,H2)
    end
    m = length(P)
    n = length(Q)
    K = exp.(-(1/eps).*c)

    function F(v)
        P.*((K*v).^(-1))
    end
    function G(u)
        Q.*((K'*u).^(-1))
    end
    #Now solve the equations u=F(v), v=G(u) by fixed-point iteration.
    u = ones(m)
    e = 1
    while e > .0001
        w = deepcopy(u)
        u = F(G(u))
        e = norm(u-w)
    end
    v = G(u)
    #compute cost = sum_{i,j} c[i,j] π[i,j] = sum_{i,j} c[i,j] u[i] K[i,j] v[i]
    cost = dot(u, (K.*c)*v)
end

#Here is a version of the Sinkhorn algorithm which takes advantage of a cancellation step.
#This can be applied iff:
#(a) The cost function satisfies c(x,x)=0 and c(x,y)<=c(x,z)+c(z,x) for all x,y,z. For example if c is a metric.
#(b) The two histograms have the same edge sets.
#Note: this function is also written ONLY for the special case that H1 and H2 are formed from samples of the same size. This is due to Julia's preference for storing Histogram weights in Int64 format.

function csolve_sinkhorn(H1::Histogram, H2::Histogram, eps::Real=0.125)
    if H1 == H2
        0.0
    else
        @assert H1.edges == H2.edges
        @assert sum(H1.weights) == sum(H2.weights)
        mydiff = H1.weights - H2.weights
        histplus = deepcopy(H1)
        histminus = deepcopy(H1)
        histplus.weights = (mydiff.>=0).*mydiff
        histminus.weights = -(mydiff.<=0).*mydiff
        lambda = sum(histplus.weights) / sum(H1.weights)
        H1 = nothing
        H2 = nothing
        cost = lambda*solve_sinkhorn(histplus, histminus, eps)
    end
end

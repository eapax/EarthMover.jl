using Hungarian
using StatsBase
include("cost_functions.jl")
include("histo_unwrap_monge.jl")

#Here is a cost-function-agnostic version of the Monge algorithm
function solve_monge(H1::Histogram, H2::Histogram)
    @assert sum(H1.weights)==sum(H2.weights)
    @assert length(H1.edges)==length(H2.edges)
    @assert length(H1.edges)<=3
    if length(H1.edges)==1
        costmatrix = cost_matrix1d(H1,H2)
    end
    if length(H1.edges)==2
        costmatrix = cost_matrix2d(H1,H2)
    end
    if length(H1.edges)==3
        costmatrix = cost_matrix3d(H1,H2)
    end
    M = sum(H1.weights)
    cost = (1/M)*hungarian(costmatrix)[2]
end

#Here is a version of the Monge algorithm which takes advantage of a cancellation step.
#This can be applied iff:
#(a) The cost function satisfies c(x,x)=0 and c(x,y)<=c(x,z)+c(z,x) for all x,y,z. For example if c is a metric.
#(b) The two histograms have the same edge sets.
function csolve_monge(H1::Histogram, H2::Histogram)
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
        cost = lambda*solve_monge(histplus, histminus)
    end
end

#Here is a version of the Hungarian Algorithm which takes as arguments some raw, unbinned data {x_i} & {y_i} and computes the Wasserstein distance between the corresponding discrete measures Sum_i delta_{x_i}

function psolve_monge(X::Array{<:Union{Tuple, Array}, 1}, Y::Array{<:Union{Tuple, Array}, 1})
    @assert length(X)==length(Y)
    M = length(X)
    c = mycost.(X, permutedims(Y))
    cost = (1/M)*hungarian(c)[2]
end

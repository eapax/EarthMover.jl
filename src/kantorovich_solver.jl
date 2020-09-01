#We will make use of the linear optimization solver provided within the scipy.optimize package.
using PyCall
spop = pyimport("scipy.optimize")
using StatsBase

include("cost_functions.jl")
include("histo_unwrap_kantorovich.jl")

#We will first write a solver which works for any mycost function.
#Note that this function automatically normalizes the histograms, if they are not already normalized

function solve_kantorovich(H1::Histogram, H2::Histogram)
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
    # Reshape the cost matrix into a vector. we will work with row-major convention, but Julia uses column-major by default. Hence the transpose.
    c = vec(c')
    # Define our constraint matrix A such that the constraints
    # \sum_j π[i,j] = P[i], \sum_i π[i,j] = Q[j]
    # can be put in canonical form \sum_r A[r,k] x[r] = b[k]
    # To this end, we will introduce the row-major form function
    function row_major(i::Int,j::Int)
        n*(i-1)+j
    end
    A = zeros(m+n,m*n)
    b = zeros(m+n)
    for k in 1:m
        for j in 1:n
            A[k,row_major(k,j)] = 1
        end
        b[k] = P[k]
    end
    for k in m+1:m+n
        for i in 1:m
            A[k,row_major(i,k-m)] = 1
        end
        b[k] = Q[k-m]
    end
    P = nothing
    Q = nothing
    #Run the linear optimization
    optimum = spop.linprog(c=c, A_eq = A, b_eq = b)
    cost = optimum["fun"]
end

#Next is a solver which takes advantage of an additional cancellation step, first taking the difference between histograms and then evaluating the
#This can be applied iff:
#(a) The cost function satisfies c(x,x)=0 and c(x,y)<=c(x,z)+c(z,x) for all x,y,z. For example if c is a metric.
#(b) The two histograms have the same edge sets.
#Note: this function is also written ONLY for the special case that H1 and H2 are formed from samples of the same size. This is due to Julia's preference for storing Histogram weights in Int64 format.

function csolve_kantorovich(H1::Histogram, H2::Histogram)
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
    cost = lambda*solve_kantorovich(histplus,histminus)
end

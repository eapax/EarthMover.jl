#The functions below recast a pair of histograms into a pair of probability vectors and a cost matrix.
#At present there are 3 distinct functions, for 1D, 2D and 3D histograms respectively. At some point I would like to change this to be dimension agnostic...

include("cost_functions.jl")

function make_probcost_vectors1d(H1::Histogram, H2::Histogram)
    @assert length(H1.edges)==length(H2.edges)==1
    edges1 = H1.edges[1]
    edges2 = H2.edges[1]

    #Create indexes of those weights which are non-zero in each histogram
    W1 = H1.weights
    W2 = H2.weights
    P = Array{Real,1}()
    Q = Array{Real,1}()

    A1 = length(W1)
    edgevect1 = Array{Real,1}()
    for i in 1:A1
        if W1[i]!=0
            append!(P, W1[i])
            append!(edgevect1, edges1[i])
        end
    end

    A2 = length(W2)
    edgevect2 = Array{Real,1}()
    for i in 1:A2
        if W2[i]!=0
            append!(Q, W2[i])
            append!(edgevect2, edges2[i])
        end
    end

    #Ensure P & Q are normalized to give probability vectors.
    P = P / sum(W1)
    Q = Q / sum(W2)
    c = mycost.(edgevect1, permutedims(edgevect2))
    return (P,Q,c)
end

function make_probcost_vectors2d(H1::Histogram, H2::Histogram)
    @assert length(H1.edges)==length(H2.edges)==2
    edges1 = H1.edges
    edges2 = H2.edges


    #Create indexes of those weights which are non-zero in each histogram
    W1 = H1.weights
    W2 = H2.weights
    P = []
    Q = []

    A1 = size(W1)[1]
    B1 = size(W1)[2]
    edgevect1 = Array{Array{Real,1},1}()
    for i in 1:A1
        for j in 1:B1
            if W1[i,j]!=0
                append!(P, W1[i,j])
                push!(edgevect1, [edges1[1][i],edges1[2][j]])
            end
        end
    end

    A2 = size(W2)[1]
    B2 = size(W2)[2]
    edgevect2 = Array{Array{Real,1},1}()
    for i in 1:A2
        for j in 1:B2
            if W2[i,j]!=0
                append!(Q, W2[i,j])
                push!(edgevect2, [edges2[1][i],edges2[2][j]])
            end
        end
    end

    #Ensure P & Q are normalized to give probability vectors.
    P = P / sum(W1)
    Q = Q / sum(W2)
    c = mycost.(edgevect1, permutedims(edgevect2))
    return (P,Q,c)
end

function make_probcost_vectors3d(H1::Histogram, H2::Histogram)
    @assert length(H1.edges)==length(H2.edges)==3
    edges1 = H1.edges
    edges2 = H2.edges

    #Create indexes of those weights which are non-zero in each histogram
    W1 = H1.weights
    W2 = H2.weights
    P = []
    Q = []

    A1 = size(W1)[1]
    B1 = size(W1)[2]
    C1 = size(W1)[3]
    edgevect1 = Array{Array{Real,1},1}()
    for i in 1:A1
        for j in 1:B1
            for k in 1:C1
                if W1[i,j,k]!=0
                    append!(P, W1[i,j,k])
                    push!(edgevect1, [edges1[1][i],edges1[2][j],edges1[3][k]])
                end
            end
        end
    end

    A2 = size(W2)[1]
    B2 = size(W2)[2]
    C2 = size(W2)[3]
    edgevect2 = Array{Array{Real,1},1}()
    for i in 1:A2
        for j in 1:B2
            for k in 1:C2
                if W2[i,j,k]!=0
                    append!(Q, W2[i,j,k])
                    push!(edgevect2, [edges2[1][i],edges2[2][j],edges2[3][k]])
                end
            end
        end
    end

    #Ensure P & Q are normalized to give probability vectors.
    P = P / sum(W1)
    Q = Q / sum(W2)
    c = mycost.(edgevect1, permutedims(edgevect2))
    return (P,Q,c)
end

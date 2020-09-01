function wcount1d(H::Histogram)
    wcount = []
    weights = H.weights
    edges = H.edges[1]
    M = size(weights)[1]
    for i in 1:M
        adj = [ edges[i] for j in 1:weights[i] ]
        wcount = vcat(wcount, adj)
    end
    return wcount
end

function wcount2d(H::Histogram)
    wcount = []
    M = size(H.weights)[1]
    N = size(H.weights)[2]
    @simd for j in 1:N
        @simd for i in 1:M
            @inbounds weight = H.weights[i,j]
            @inbounds bin_coords = (H.edges[1][i],H.edges[2][j])
            @simd for l in 1:weight
                wcount = vcat(wcount, bin_coords)
            end
        end
    end
    return wcount
end

function wcount3d(H::Histogram)
    wcount = []
    M = size(H.weights)[1]
    N = size(H.weights)[2]
    P = size(H.weights)[3]
    @simd for k in 1:P
        @simd for j in 1:N
            @simd for i in 1:M
                @inbounds weight = H.weights[i,j,k]
                @inbounds bin_coords = (H.edges[1][i],H.edges[2][j],H.edges[3][k])
                @simd for l in 1:weight
                    wcount = vcat(wcount, bin_coords)
                end
            end
        end
    end
    return wcount
end

function cost_matrix1d(H1::Histogram, H2::Histogram)
    x = wcount1d(H1)
    y = permutedims(wcount1d(H2))
    return mycost.(x,y)
end

function cost_matrix2d(H1::Histogram, H2::Histogram)
    x = wcount2d(H1)
    y = permutedims(wcount2d(H2))
    return mycost.(x,y)
end

function cost_matrix3d(H1::Histogram, H2::Histogram)
    x = wcount3d(H1)
    y = permutedims(wcount3d(H2))
    return mycost.(x,y)
end

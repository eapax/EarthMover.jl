using LinearAlgebra

function mycost(x::Union{Array{<:Real}, Real, Tuple},y::Union{Array{<:Real}, Real, Tuple})
    @assert length(x)==length(y)
    answer = norm(x.-y)
end

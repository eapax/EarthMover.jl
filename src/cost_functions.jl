using LinearAlgebra

function mycost(x,y)
    @assert length(x)==length(y)
    answer = norm(x.-y)
end

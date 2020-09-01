module EarthMover

export solve_kantorovich, solve_monge, solve_sinkhorn

include("kantorovich_solver.jl")
include("monge_solver.jl")
include("sinkhorn_solver.jl")

end # module

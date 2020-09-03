using EarthMover, StatsBase, Test

#The following 3d example can be solved by hand, giving Wass1(H1,H2)=0.75
@testset "First 3D exactly computable example" begin
    data1 = ([1,2,4,4],[2,2,4,1],[1,1,1,2])
    data2 = ([1,3,2,4],[2,4,2,2],[1,1,2,2])
    edges = (1:5,1:5,1:3)
    epsilon = 0.01
    H1 = fit(Histogram, data1, edges)
    H2 = fit(Histogram, data2, edges)
    x = solve_kantorovich(H1,H2)
    #First test positive definiteness
    @test solve_kantorovich(H1,H1)≈0.0 atol=0.0001 
    @test csolve_kantorovich(H1,H1)==0.0
    @test solve_monge(H1,H1)==0.0
    @test csolve_monge(H1,H1)==0.0
    @test solve_sinkhorn(H1,H1)≈0.0 atol=0.01
    @test csolve_sinkhorn(H1,H1)==0.0
    #Now test for exact values
    @test x ≈ 0.75 atol=0.001
    @test x ≈ csolve_kantorovich(H1,H2)  atol=0.001
    @test x ≈ solve_monge(H1, H2)  atol=0.001
    @test solve_monge(H1, H2) ≈ csolve_monge(H1, H2)  atol=0.001
    @test x ≈ solve_sinkhorn(H1, H2, epsilon)  atol=0.01
    @test solve_sinkhorn(H1, H2, epsilon) ≈ csolve_sinkhorn(H1, H2, epsilon)  atol=0.001
end

#The following 3d example can also be solved by hand, giving Wass1(H1,H2)=2/3
@testset "Second 3D exactly computable example" begin
    data1 = ([.5,1.5,1.5],[.5,1.5,2.5],[.5,.5,.5])
    data2 = ([.5,2.5,.5],[.5,1.5,2.5],[.5,.5,.5])
    edges = (0:3,0:3,0:3)
    H1 = fit(Histogram, data1, edges)
    H2 = fit(Histogram, data2, edges)
    epsilon = 0.01
    x = solve_kantorovich(H1,H2)
    @test x ≈ 2/3 atol=0.001
    @test x ≈ csolve_kantorovich(H1,H2)  atol=0.001
    @test x ≈ solve_monge(H1, H2)  atol=0.001
    @test solve_monge(H1, H2) ≈ csolve_monge(H1, H2)  atol=0.001
    @test x ≈ solve_sinkhorn(H1, H2, epsilon)  atol=0.01
    @test solve_sinkhorn(H1, H2, epsilon) ≈ csolve_sinkhorn(H1, H2, epsilon) atol=0.001
end

#Actually, the above example was really 2D, so lets also check the 2D solver with it
@testset "2D exactly computable example" begin
    data1 = ([.5,1.5,1.5],[.5,1.5,2.5])
    data2 = ([.5,2.5,.5],[.5,1.5,2.5])
    edges = (0:3,0:3)
    H1 = fit(Histogram, data1, edges)
    H2 = fit(Histogram, data2, edges)
    epsilon = 0.01
    x = solve_kantorovich(H1,H2)
    @test x ≈ 2/3 atol=0.001
    @test x ≈ csolve_kantorovich(H1,H2)  atol=0.001
    @test x ≈ solve_monge(H1, H2)  atol=0.001
    @test solve_monge(H1, H2) ≈ csolve_monge(H1, H2)  atol=0.001
    @test x ≈ solve_sinkhorn(H1, H2, epsilon)  atol=0.01
    @test solve_sinkhorn(H1, H2, epsilon) ≈ csolve_sinkhorn(H1, H2, epsilon) atol=0.001
end

#Now for a 1D example. It can be checked by hand that Wass1(H1,H2)=0.12. This one is important as it is the first non-square problem for Kantorovich (i.e m!=n)
@testset "1D exactly computable example" begin
    data1 = [.1,.1,.1,.3,.3,.3,.3,.5,.5,.5]
    data2 = [.1,.1,.3,.3,.3,.3,.5,.7,.7,.9]
    edges = 0:.2:1
    epsilon = 0.01
    H1 = fit(Histogram, data1, edges)
    H2 = fit(Histogram, data2, edges)
    x = solve_kantorovich(H1,H2)
    @test x ≈ 0.12 atol=0.001
    @test x ≈ csolve_kantorovich(H1,H2)  atol=0.001
    @test x ≈ solve_monge(H1, H2)  atol=0.001
    @test solve_monge(H1, H2) ≈ csolve_monge(H1, H2)  atol=0.001
    @test x ≈ solve_sinkhorn(H1, H2, epsilon)  atol=0.01
    @test solve_sinkhorn(H1, H2, epsilon) ≈ csolve_sinkhorn(H1, H2, epsilon) atol=0.001
end

#Now some randomized examples for good measure!

#1D
@testset "1D randomized example" begin
    data1 = rand(100)
    data2 = rand(100)
    edges = 0:.2:1
    H1 = fit(Histogram, data1, edges)
    H2 = fit(Histogram, data2, edges)
    epsilon = 0.05
    x = solve_kantorovich(H1,H2)
    @test x ≈ csolve_kantorovich(H1,H2)  atol=0.01 rtol=0.01
    @test x ≈ solve_monge(H1, H2)  atol=0.01 rtol=0.01
    @test solve_monge(H1, H2) ≈ csolve_monge(H1, H2)  atol=0.01 rtol=0.01
    @test x ≈ solve_sinkhorn(H1, H2, epsilon)  atol=0.2 rtol=0.1
    @test solve_sinkhorn(H1, H2, epsilon) ≈ csolve_sinkhorn(H1, H2, epsilon) atol=0.2 rtol=0.1
end

#2D
@testset "2D randomized example" begin
    data1 = (rand(100), rand(100))
    data2 = (rand(100), rand(100))
    edges = (0:.2:1, 0:.2:1)
    H1 = fit(Histogram, data1, edges)
    H2 = fit(Histogram, data2, edges)
    epsilon = 0.05
    x = solve_kantorovich(H1,H2)
    @test x ≈ csolve_kantorovich(H1,H2)  atol=0.01 rtol=0.01
    @test x ≈ solve_monge(H1, H2)  atol=0.01 rtol=0.01
    @test solve_monge(H1, H2) ≈ csolve_monge(H1, H2)  atol=0.01 rtol=0.01
    @test x ≈ solve_sinkhorn(H1, H2, epsilon)  atol=0.2 rtol=0.1
    @test solve_sinkhorn(H1, H2, epsilon) ≈ csolve_sinkhorn(H1, H2, epsilon) atol=0.2 rtol=0.1
end

#3D
@testset "3D randomized example" begin
    data1 = (rand(100), rand(100), rand(100))
    data2 = (rand(100), rand(100), rand(100))
    edges = (0:.2:1, 0:.2:1, 0:.2:1)
    H1 = fit(Histogram, data1, edges)
    H2 = fit(Histogram, data2, edges)
    epsilon = 0.05
    x = solve_kantorovich(H1,H2)
    @test x ≈ csolve_kantorovich(H1,H2)  atol=0.01 rtol=0.01
    @test x ≈ solve_monge(H1, H2)  atol=0.01 rtol=0.01
    @test solve_monge(H1, H2) ≈ csolve_monge(H1, H2)  atol=0.01 rtol=0.01
    @test x ≈ solve_sinkhorn(H1, H2, epsilon)  atol=0.2 rtol=0.1
    @test solve_sinkhorn(H1, H2, epsilon) ≈ csolve_sinkhorn(H1, H2, epsilon) atol=0.2 rtol=0.1
end

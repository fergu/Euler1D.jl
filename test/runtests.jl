using Test
using Euler1D

# Test the callback functions
@testset "Callbacks" begin
    @testset "Cycle Callbacks" begin
        include("CycleCallbacks.jl")
    end
    @testset "Time Callbacks" begin
        include("TimeCallbacks.jl")
    end
    @testset "Time Delta Callbacks" begin
        include("TimeDeltaCallbacks.jl")
    end
end

# Test Sod's shock tube
@testset "Sod Shock Tube" begin
    include("SodShockTube.jl")
end

using Test
using Euler1D

# Test the callback functions
@testset "Callbacks" begin
    include("CycleCallbacks.jl")
end

# Test Sod's shock tube
@testset "Sod Shock Tube" begin
    include("SodShockTube.jl")
end

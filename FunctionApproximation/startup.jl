using Pkg; Pkg.activate(".")
using Distributed
using Plots
using Distributions
using GLMakie 
@everywhere begin
	using Pkg; Pkg.activate(".")
	Pkg.instantiate(); Pkg.precompile()
end
@everywhere begin
	using ProgressMeter
end
@everywhere begin
	include("function.jl")
end




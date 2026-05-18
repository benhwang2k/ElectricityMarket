using Roots


include("startup.jl")
using JSON
using Random



"""
The idea here is to model a virtual power plant (vpp) that 
consists of an energy storage device (battery) that is 
connected to both the electrical grid and a solar cell.

The battery is used to serve a random net energy signal
that is fluctuates randomly from demand (when grid needs energy)
to supply (when solar is higher than demand). 

Idea: multiple batteries may be in the vpp with different
distributions of data. We can model these independently 
because maybe they will have different optimal control inputs
so we can have them as different dimensions in the state.
 IE
[x1, x2] = [soc_1, soc_2]  are the different states of charge 
of the different batteries that recieve different random inputs of net demand.

this demand could be from the same distribution.

"""

M = 20
grid = [[i/M for i in 0:M] for _ in 1:1]
cov = 0.2#0.2 #0.1
alpha = 0.8 


# Transition dynmics
function reflect(x)
	count = 0
	while (x < 0) || (x > 1)
		if x < 0
			x = -1 * x
		end
		if x > 1
			x = 2 - x
		end
		count += 1 
	end
	if count > 1
		println("reflection count = $count ")
	end
	return x
end


demand_data = JSON.parsefile("demand_data.json")
solar_data = JSON.parsefile("solar_data.json")

data_d =  [parse(Int, d["value"]) for d in demand_data if d["value"] !== nothing] 
data_s =  [parse(Int, d["value"]) for d in solar_data if d["value"] !== nothing]
n_data_d = Int(floor(length(data_d) / 24))
data_d = vec(sum(reshape(data_d[1:(24*n_data_d)], 24, :), dims=1))
n_data_s = Int(floor(length(data_s) / 24))
data_s = vec(sum(reshape(data_s[1:(24*n_data_s)], 24, :), dims=1))


data_d = data_d .- minimum(data_d)
data_d = data_d .* (1.0/maximum(data_d))
mean_d = sum(data_d) / length(data_d)
data_s = data_s .- minimum(data_s)
data_s = data_s .* (1.0/maximum(data_s))
mean_s = sum(data_s) / length(data_s)


# These boundary conditions dont result in a Gaussian
#
function sample_next(x)
	y = x[1] + (rand(data_s) - mean_s - rand(data_d) + mean_d)/10.0 
	y = min(y, 1.0)
	y = max(y, 0.0)
	return y
end

# Initial Condition
x = 0.5

# Run it for data
println("creating data")
data = [[x]]
D = 6 
for t in 1:(10^D)
	push!(data, [sample_next(data[end])])
end

# If you want the histogram version of the data;
do_histogram = true 
if do_histogram
	println("plotting histogram")
	datah = [d[1] for d in data]
	ph = Plots.histogram(datah, normalize= :pdf, label="Histogram")
	display(ph)
end

#=

# Here is the Theoretical Stationary Distribtution
gram = cov^2 * (1 / (1-(1-alpha)^2)) 
plotpdf = Normal(0.5, sqrt(gram))
renorm_fac = cdf(plotpdf,1.0) - cdf(plotpdf, 0.0)

function thr_pdf(x) 
	return pdf(plotpdf, x) + pdf(plotpdf, -1*x) + pdf(plotpdf, 2-x) 
end

function pdf_derivative(x)
	return (1 / sqrt(2*pi*gram))*(-2*(x-0.5)/(2*gram))*exp((-1*(x-0.5)^2)/(2*gram))
end

function pdf_derivative2(x)
	return (1 / sqrt(2*pi*gram))*((((x-0.5)^2)/(gram^2))-(1/gram))*exp((-1*(x-0.5)^2)/(2*gram))
end

function reflected_derv(x)
	return pdf_derivative(x) - pdf_derivative(-1*x) - pdf_derivative(2-x)
end

function reflected_derv2(x)
	return pdf_derivative2(x) + pdf_derivative2(-1*x) +  pdf_derivative2(2-x)
end





Plots.plot!(thr_pdf, 0.0:0.001:1.0)
display(ph)

=#
# Derivative nonsense
#=
pd = Plots.plot(reflected_derv, 0:0.01:1)
rt1 = find_zero(reflected_derv2, 0.5 - sqrt(gram))
rt2 = find_zero(reflected_derv2, 0.5 + sqrt(gram))
println("rts = $([rt1, rt2])")
Plots.plot!([rt1], [reflected_derv(rt1)], markershape=:circle)
Plots.plot!([rt2], [reflected_derv(rt2)], markershape=:circle)
println("rts = $([rt1, rt2])")
display(pd)
=#

# Here is the function for the kernel density
function kern_dens(x, y)
	# this is a gaussian reflected at the borders 0, 1
	random_var = x + alpha*(0.5 - x) + cov*Normal(0,1) 
	return pdf(random_var, y) + pdf(random_var, -1*y) + pdf(random_var, 2-y)
end

#finegrid = [i for i in 0:0.01:1]
#kern = [kern_dens(x, y) for x in finegrid, y in finegrid]
#pheat_true = Plots.heatmap(finegrid, finegrid, kern)



println("type 1 to continue. Else exit")
if true#r == "1"

	# Frankliin Approximation!

	# 1D basis for the stationary distribution
	println("1D basis approximation")
	numbasis = 2^5 + 1#2^6
	basis = @time Func.create_basis(1, numbasis)
	println("time for emp app")
	em1_c = @time Func.approximate_emp(basis, data) 
	em1 = em1_c[1]
	cem1 = em1_c[2]

	nem1 = Func.normalize_1(em1) 

	# plot the approximation of the stationary
	pn = Plots.plot!(x -> Func.eval_f(nem1, [x]), 0:0.01:1, label=L"\hat{P}_\pi")
	Plots.xlabel!("State")
	Plots.ylabel!("Probability Density")
	display(pn)
end
println("type 1 to continue. Else exit")
if false#r == "1"
	# Use the stationary approximation for the kernel approximation

	# construct the paired dataset
	data2 = [[data[i][1], data[i+1][1]] for i in 1:length(data)-1]
	println("time for 2d basis")
	basis2 = @time Func.create_basis(2, 8)#16)

	println("time for kern app")
	
	em2_c = @time Func.approximate_kern(basis2, data2, em1)
	em2 = em2_c[1]
	cem2 = em2_c[2]

	em2 = Func.normalize_1(em2)
	
	z = Func.plot(em2)
	Makie.surface(z...)
	println(minimum(z[3]))

end



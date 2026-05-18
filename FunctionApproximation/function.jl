module Func 
import Base:-,+,*
using DynamicPolynomials
using Distributed
using Pkg
using ProgressMeter
using MultiQuad
using QuadGK
using Serialization
export eval_fnc, fnc, MultiLinearFunction, PWMLinearFunction 

abstract type fnc end


function getCombs(x1::Vector{Float64}, x2::Vector{Float64})
	# input: two corners
	# output all corners of the box defined by the corners
	if size(x1) != size(x2)
		error("inputs just have same dim")
	end
	if size(x1)[1] == 1 && size(x2)[1] == 1
		return [Float64[x1[1]], Float64[x2[1]]]
	end
	combs = getCombs(x1[2:end], x2[2:end])
	newcombs = Vector{Float64}[] 
	for c in combs
		nc1 = copy(c)
		nc2 = copy(c)
		pushfirst!(nc1, x1[1])
		pushfirst!(nc2, x2[1])
		push!(newcombs, nc1)
		push!(newcombs, nc2)
	end
	return newcombs
end

struct MultiLinearFunction <: fnc
	x1::Vector{Float64}
	x2::Vector{Float64}
	y::Dict{Vector{Float64}, Float64}
	function MultiLinearFunction(x1::Vector{Float64}, x2::Vector{Float64}, y::Dict{Vector{Float64}, Float64})
		if x1 >= x2
			error("x1 must be less than x2")
		end
		# iterate through and make sure that the dict has all values	
		new(copy(x1), copy(x2), deepcopy(y))
	end
end


struct PWMLinearFunction <: fnc 
	pieces::Dict{Vector{Float64}, MultiLinearFunction}
end

function getGrid(f::PWMLinearFunction)
	dim = length(first(keys(f.pieces)))
	xs = [unique(sort([pt[i] for pt in keys(f.pieces)])) for i in 1:dim]
	return xs
end

function Base.:*(c::Float64, f::MultiLinearFunction)
	newd = Dict{typeof(first(keys(f.y))), typeof(f.y[first(keys(f.y))])}()
	sizehint!(newd, length(f.y))             # avoid rehashing

	@inbounds for k in keys(f.y)
		newd[k] = c * f.y[k]
	end

	return MultiLinearFunction(f.x1, f.x2, newd)
end


function Base.:+(f::MultiLinearFunction, g::MultiLinearFunction)
	newd = Dict{typeof(first(keys(f.y))), typeof(f.y[first(keys(f.y))])}()
    sizehint!(newd, length(f.y))
 
    @inbounds for k in keys(f.y)
        newd[k] = f.y[k] + g.y[k]
    end

    return MultiLinearFunction(f.x1, f.x2, newd)
end

function abs(f::MultiLinearFunction)
	newd = Dict{typeof(first(keys(f.y))), typeof(f.y[first(keys(f.y))])}()
    sizehint!(newd, length(f.y))
 
    @inbounds for k in keys(f.y)
	    newd[k] = Base.abs(f.y[k])
    end

    return MultiLinearFunction(f.x1, f.x2, newd)
end

Base.:-(f::MultiLinearFunction, g::MultiLinearFunction) = f + (-1.0 * g)


function Base.:*(c::Float64, f::PWMLinearFunction)
	newd = Dict{typeof(first(keys(f.pieces))), typeof(f.pieces[first(keys(f.pieces))])}()
	sizehint!(newd, length(f.pieces))
	
	for k in keys(f.pieces)
		newd[k] = c * f.pieces[k]
	end

	return PWMLinearFunction(newd)
end

function abs(f::PWMLinearFunction)
	newd = Dict{typeof(first(keys(f.pieces))), typeof(f.pieces[first(keys(f.pieces))])}()
	sizehint!(newd, length(f.pieces))

	for k in keys(f.pieces)
		newd[k] = abs(f.pieces[k])
	end

	return PWMLinearFunction(newd)
end

function Base.:+(f::PWMLinearFunction, g::PWMLinearFunction)
	newd = Dict{typeof(first(keys(f.pieces))), typeof(f.pieces[first(keys(f.pieces))])}()
	sizehint!(newd, length(f.pieces))

	for k in keys(f.pieces)
		newd[k] = f.pieces[k] + g.pieces[k]
	end

	return PWMLinearFunction(newd)
end

function Base.:+(f::PWMLinearFunction, g::Function)
	newd = Dict{typeof(first(keys(f.pieces))), typeof(f.pieces[first(keys(f.pieces))])}()
	sizehint!(newd, length(f.pieces))
	for k in keys(f.pieces)
		newd[k] = f.pieces[k] + g(k...)
	end
	return PWMLinearFunction(newd)
end

Base.:+(g::Function, f::PWMLinearFunction) = f + g
Base.:-(g::PWMLinearFunction, f::Function) = (-1.0*f) + g
Base.:-(g::Function, f::PWMLinearFunction) = f + (-1.0*g)

Base.:-(f::PWMLinearFunction, g::PWMLinearFunction) = f + (-1.0 * g)

# Constructor for  PWLinearFunction
function PWMLinearFunction(fs::Vector{MultiLinearFunction})
	d = Dict{Vector{Float64}, MultiLinearFunction}()
	for f in fs 
		d[f.x1] = f
	end
	return PWLinearFunction(d)
end



function plot(f::PWMLinearFunction)
    println("plotting…")

    pts = collect(keys(f.pieces))
    dim = length(pts[1])   # number of coordinates in each point

    if dim == 1
        # -------------------------
        # 1D PLOT
        # -------------------------
	x2 = f.pieces[maximum([c for c in keys(f.pieces)])].x2
        xs = sort(unique([p[1] for p in pts]))
	append!(xs, x2[1])
	ys = [eval_f(f, [x]) for x in xs]

        return (xs, ys)

    elseif dim == 2
        # -------------------------
        # 2D SURFACE PLOT
        # -------------------------
	grid = getGrid(f)
	x2 = f.pieces[maximum([c for c in keys(f.pieces)])].x2
        xs = sort(unique([p[1] for p in pts]))
	append!(xs, x2[1])
        ys = sort(unique([p[2] for p in pts]))
	append!(ys, x2[2])

        #Z = [f.pieces[[x, y]].y[[x, y]] for x in xs, y in ys]
	Z = [eval_f(f,[x,y]) for x in xs, y in ys]

        return (xs, ys, Z)

    else
        # -------------------------
        # ND SURFACE PLOT
        # -------------------------
	grid = getGrid(f)
	x2 = f.pieces[maximum([c for c in keys(f.pieces)])].x2
	i = 1
	for g in grid
		append!(g, x2[i])
		i += 1
	end
	Z = Array{Float64}(undef, [length(g) for g in grid]...) 
	for i in CartesianIndices(Z)
		pt = [grid[j][i[j]] for j in 1:length(grid)]
		Z[i] = eval_f(f, pt) 
	end
        return (grid, Z)
        #error("Plotting for dimension $dim is not supported.")
    end
end

function supnorm(f::MultiLinearFunction)
	return max([Base.abs(f.y[k]) for k in keys(f.y)])
end

function supnorm(f::PWMLinearFunction)
	return max([supnorm(f.pieces[c]) for c in keys(f.pieces)])
end

function PWMLinearFunction(grid::Vector{Vector{Float64}}, f)
	# grid is a vector of discretized axes [[0,0.5,1],[0,0.5,1]]
	smallgrid = [v[1:end-1] for v in grid]
	dim = length(grid)
	d = Dict{Vector{Float64}, MultiLinearFunction}()
	for idx in Iterators.product((eachindex(x) for x in smallgrid)...)
		x1 = [ grid[k][idx[k]] for k in 1:length(grid) ]
		x2 = [ grid[k][idx[k] + 1] for k in 1:length(grid) ]
		dc = Dict([c => f(c) for c in getCombs(x1, x2)])
		fnew = Func.MultiLinearFunction(x1, x2, dc) 
		d[x1] = fnew
	end
	return PWMLinearFunction(d)
end

	



function eval_fnc(f::MultiLinearFunction, x::Vector{Float64})
	if size(x,1) != size(f.x1, 1) 
		error("point must be of same dimension as function domain") 
	end 
	if !(f.x1 <= x <= f.x2) 
		error("point must be between limits of function. x = $x, x1 = $x1, x2 = $x2") 
	end 
	val = 0 
	# The next code parameterizes every convex combintation of the values at the corners in terms of the position in the cell. 
	for corner in keys(f.y)
		prodresult = 1 
		for n in 1:size(f.x1,1)
			t = (x[n] .- f.x1[n])/(f.x2[n] .- f.x1[n]) 
			prodresult *= corner[n] == f.x1[n] ? 1-t : t
		end 
		val += f.y[corner] * (prodresult) 
	end 
	return val 
end

function eval_f(f::PWMLinearFunction, x::Vector{Float64}) 
	grid = getGrid(f)
	dim = length(x)
	ll = Vector{Float64}(undef, dim)
	for d in 1:dim
		gd = grid[d]
		i = searchsortedlast(gd, x[d])
		# clamp index to valid range
		i = max(1, min(i, length(gd)))
		ll[d] = gd[i]
	end
	val = -1
	try 	
		val = eval_fnc(f.pieces[ll], x) 
	catch e
		println("ERROR evaluating function \n at point $x \n \n \n")
		println(f)
		println(x)
		println(ll)
		println(e)
	end
	return val
end

function weight(c,x)
	p = 1
	for i in eachindex(c)
		p *= (c[i] == 0 ? (1 - x[i]) : x[i])
	end
	return p
end

#=
function refine(f::PWMLinearFunction)
	grid = getGrid(f)
	n = length(grid[1])
	newg = Vector
=#

function integrate(f::MultiLinearFunction, x1, x2, x)
	# This function integrates f over the box [x1, x2]
	# If x1 < f.x1 or x2 > f.x2, then the integral is take only over the domain of f
	# If x1 > f.x2 or x2 < f.x1 returns 0
	dimension = length(first(keys(f.y)))
	m = 0
	if x1 > f.x2
		return 0
	end
	if x2 < f.x1
		return 0
	end
	if x1 < f.x1
		x1 = f.x1
	end
	if x2 > f.x2
		x2 = f.x2
	end
	width = x2 - x1
	corners = getCombs(zeros(Float64, dimension), ones(Float64, dimension))
	corners = [[round(Int64, y) for y in x] for x in corners]
	for c1 in corners
		p = weight(c1, x)
		cval1 = eval_fnc(f, x1 + width .* c1) 
		for i in 1:dimension
			p = antidifferentiate(p, x[i]) 
		end
		m += subs(p, x => width) == 0.0 ? 0.0 : (cval1 * coefficients(subs(p, x => width))[1])
	end
	return m
end

function integrate(f::PWMLinearFunction, x1, x2)
	dimension = length(first(keys(f.pieces)))
	@polyvar x[1:dimension]
	total = 0
	keylist = collect(keys(f.pieces))
	l = length(keylist)
	n = nworkers()
	split = [round(Int, (l*i/n)) for i in 0:n]
	chunks = [keylist[split[k] + 1 : split[k+1]] for k in 1:n] 
	if nworkers() > 1
		totals = pmap(chunks) do c
			total = 0
			for k in c 
				total += integrate(f.pieces[k], x1, x2, x)
			end
			total 
		end
	else
		totals = [0.0]
		for k in keylist
			totals[1] += integrate(f.pieces[k], x1, x2, x)
		end
	end
	return sum(totals) 
end

function integrate_upto(f::PWMLinearFunction, x)
	return integrate(f, sort(collect(keys(f.pieces)))[1], x)
end

function integrate_all(f::PWMLinearFunction)
	c2 = maximum([c for c in keys(f.pieces)])
	x2 = f.pieces[c2].x2
	return integrate_upto(f, x2)
end


function slice(f::PWMLinearFunction, x::Float64; d=1)
	# This function slices off the d dimension at value x
	# and returns the lower dimensional PWMLinearFunction 
	# evaluated at [x, .... ]
	grid = getGrid(f)
	if !(1 <= d <= length(grid))
		error("d must be less than dimension of the function")
	end
	if length(grid) == 1
		return eval_f(f, [x])
	end
	grid = grid[[1:d-1; d+1:end]]
	return PWMLinearFunction(grid, y -> eval_f(f, insert!(y,d,x))) 
end

function my_find_zero(f, x, lower=0.0, upper=1.0, tol=0.001)
	if !(lower < x < upper)
		error("x must be inbetween lower and upper")
	end
	while upper - lower > tol
		if f(x) > 0
			upper = x
		else
			lower = x
		end
		x = (lower + upper) /2.0
	end
	return x
end


function invert_cdf_1D(p::PWMLinearFunction, u::Float64)
	firstk = first(sort(collect(keys(p.pieces))))
	midk = sort(collect(keys(p.pieces)))[ceil(Int, length(collect(keys(p.pieces)))/2)] 
	lastk = p.pieces[sort(collect(keys(p.pieces)))[end]].x2
	integ = integrate_upto(p, lastk)
	if abs(integ - 1) > 0.01
		println("p doesnt integrate close to 1 : $integ check that it is a proper density")
	end
	if u > integ 
		return lastk[1]
	end
	if u < 0.001
		return firstk[1]
	end
	dim = length(first(keys(p.pieces))) 
	if length(first(keys(p.pieces))) == 1
		return my_find_zero(x -> integrate_upto(p, Float64[x]) - u, midk[1]) 
	else
		return my_find_zero(x -> integrate_upto(p, [Float64[x]; lastk[2:end]]) - u, midk[1]) 
	end
end

function sample_pdf(p::PWMLinearFunction)
	dimension = length(first(keys(p.pieces)))
	x = [0.0 for _ in 1:dimension]
	f = p
	for d in 1:dimension
		u = rand()
		#println("dim $d : u = $u")
		x[d] = invert_cdf_1D(f, u)
		if d < dimension
			#println("dim $d : x[d] = $(x[d])")
			f = slice(f, x[d])
			lastk = f.pieces[sort(collect(keys(f.pieces)))[end]].x2
			intf = Func.integrate_upto(f, lastk)
			f = (1/intf)*f
		end
	end
	#print(x)
	return x
end




function inner_box_x(f::MultiLinearFunction, g::MultiLinearFunction, x)
	dimension = length(first(keys(f.y)))
	m = 0 
	width = f.x2 - f.x1
	corners = getCombs(zeros(Float64, dimension), ones(Float64, dimension))
	corners = [[round(Int64, y) for y in x] for x in corners]
	p = 0.0
	for c1 in corners
		W1 = weight(c1, x)
		cval1 = f.y[f.x1 + width .* c1]
		for c2 in corners
			W2 = weight(c2 ,x)
			#println("cval1 = $cval1, cval2 = $cval2, ")#W1 = $W1, W2 = $W2")
			cval2 = g.y[g.x1 + width .* c2]
			p += cval1*cval2*W1 * W2
			#=
			for i in 1:dimension
				p = antidifferentiate(p, x[i]) 
			end
			#println("intgral form = $p")
			println("width = $width, cval1 = $cval1, cval2 = $cval2, endval = $((cval1 * cval2 * coefficients(subs(p, x => width))[1]))")
			m += (cval1 * cval2 * coefficients(subs(p, x => width))[1])
			=#
		end
	end
	expansion = prod(width)
	for i in 1:dimension
		p = antidifferentiate(p, x[i]) 
	end
	covec = coefficients(subs(p, x => [1.0 for _ in 1:dimension]))
	if length(covec) == 0
		if p(x => [1.0 for _ in 1:dimension]) == 0.0
			return 0.0
		end
	end
	m = (coefficients(subs(p, x => [1.0 for _ in 1:dimension]))[1])*(expansion)
	return m
end

function inner_box(f::MultiLinearFunction, g::MultiLinearFunction)
	dimension = length(first(keys(f.y)))
	@polyvar x[1:dimension]
	return inner_box_x(f,g,x)
end


function *(f::PWMLinearFunction, g::PWMLinearFunction)
	total = 0
	dimension = length(first(keys(f.pieces[first(keys(f.pieces))].y)))
	#println("from * : dimension = $dimension")
	@polyvar x[1:dimension]
	if nworkers() > 1 && length(keys(f.pieces)) > 10 * nworkers()
		keylist = collect(keys(f.pieces))
		s = ceil(Int, length(keylist) / nworkers())
		chunks = [keylist[s*m + 1 : min(s*(m+1), length(keylist))] for m in 0:nworkers() - 1]
		result = pmap(chunks) do listi
			subtot = 0
			for k in listi
				subtot += inner_box_x(f.pieces[k], g.pieces[k], x)
			end
			return subtot
		end
		total = sum(result)
	else
		for k in keys(f.pieces)
			total += inner_box_x(f.pieces[k], g.pieces[k], x)
		end
	end
	return total
end

function *(f::PWMLinearFunction, g::Function)
	if length(first(keys(f.pieces))) > 1
		error("dimension greater than 1 not yet supported")
	end
	x1 = minimum(keys(f.pieces))
	x2 = f.pieces[maximum(keys(f.pieces))].x2
	res = QuadGK.quadgk(x -> eval_f(f, [x]) * g(x), x1[1], x2[1], maxevals=10000)
	if Base.abs(res[2]) > 0.00001	
		println("numerical integration error measurment = $(res[2])")
	end
	return round(res[1]; digits=4) 
end

function normalize_1(f::PWMLinearFunction)
	g = Func.abs(f)
	return (1.0/integrate_all(g))*g
end

function norm_sqr(corner_values, width)
	dimension = length(first(keys(corner_values)))
	m = 0 
	corners = getCombs(zeros(Float64, dimension), ones(Float64, dimension))
	corners = [[round(Int64, y) for y in x] for x in corners]
	@polyvar x[1:dimension]
	for c1 in corners
		W1 = weight(c1, x)
		cval1 = corner_values[c1]
		for c2 in corners
			W2 = weight(c2 ,x)
			cval2 = corner_values[c2]
			p = W1 * W2
			for i in 1:dimension
				p = antidifferentiate(p, x[i]) 
			end
			#m += (cval1 * cval2 * subs(p, x => width))
			m += cval1 * cval2 * p
		end
	end
	return m
end

function create_coefficient_map(dimension::Int64, width)
	m = Dict{Tuple{Vector{Int64}, Vector{Int64}}, Any}()
	corners = getCombs(zeros(Float64, dimension), ones(Float64, dimension))
	corners = [[round(Int64, y) for y in x] for x in corners]
	@polyvar x[1:dimension]
	for c1 in corners
		W1 = weight(c1, x)
		for c2 in corners
			W2 = weight(c2 ,x)
			p = W1 * W2
			for i in 1:dimension
				p = antidifferentiate(p, x[i]) 
			end



			m[c1, c2] = p 
		end
	end
	return x, m
end

function create_coefficient_map_dist(dimension::Int64)
	m = Dict{Tuple{Vector{Int64}, Vector{Int64}}, Any}()
	corners = getCombs(zeros(Float64, dimension), ones(Float64, dimension))
	corners = [[round(Int64, y) for y in x] for x in corners]
	@polyvar x[1:dimension]
	Wmap = Dict()
	status = @showprogress pmap(corners) do i
		try
			i => weight(i, x)
		catch e
			false
		end
	end
	m = Dict(status)
	n = nworkers()     # number of workers (NOT counting the main process)
	chunk_size = ceil(Int, length(corners) / n)
	chunks = [ corners[(i-1)*chunk_size+1 : min(i*chunk_size, length(corners))] for i in 1:n ]
	status = @showprogress pmap(chunks; on_error=println) do listi 
		ps = []
		num = 1
		for i in listi
			for c2 in corners
				p = m[i] * m[c2] 
				for j in 1:dimension
					p = antidifferentiate(p, x[j]) 
				end
				append!(ps, [(i, c2) => coefficients(p)]) 
			end
			num += 1
		end
		ps
	end
	return x, status  
end



struct BasisFunc
    a::Vector{Float64}
end

# Make BasisFunc callable: v(x)
function (f::BasisFunc)(x)
    tot = 1.00
    @inbounds for c in eachindex(f.a)
        ac = f.a[c]
        if ac != -1
            xc = x[c]
            tot *= (xc < ac) ? 0.00 : (xc - ac)
        end
    end
    return tot
end


function create_basis(dim, M; load_from_file=true)
	if load_from_file 
		try 
			return deserialize("basis$(dim)_$(M).jls")
		catch e
			if isa(e, SystemError)
				println("basis does not exist yet. creating ...")
			else
				error(e)
			end
		end
	end
					
	total = M^dim
	gridsize = 2^(ceil(Int, log(2,M)))
	# Preallocate
	n = fill(-1.0, dim)
	d = fill(1.0,  dim)
	index = Vector{Int}(undef, dim)
	a = Vector{Float64}(undef, dim)

	grid = [Float64[i/gridsize for i in 0:gridsize] for _ in 1:dim]
	println("grid size  = $([length(g) for g in grid]) ")
	# Preallocate basis with correct type
	basis = Vector{Func.PWMLinearFunction}(undef, 0)
	# Progress bar
	p = Progress(total; desc="Forming Basis")
	for i in 0:(total-1)
		# multi-index
		for c in 1:dim
			index[c] = div(i, M^(c-1)) % M
			a[c] = n[c] / d[c]
		end
		# update counters
		update_index = 1
		while update_index < dim && index[update_index] == M-1
			update_index += 1
			for r in 1:(update_index-1)
				n[r] = -1
				d[r] = 1
			end
		end

		if n[update_index] == -1
			n[update_index] = 0
		elseif n[update_index] + 1 == d[update_index]
			n[update_index] = 1
			d[update_index] *= 2
		else
			n[update_index] += 2
		end

		# Create basis function
		v = BasisFunc(copy(a))     # copy because a mutates!
		g = Func.PWMLinearFunction(grid, v)
		

		# Gram-Schmidt orthonormalization
		for b in basis
			g -= (b*g)*b
		end

		# normalize
		gn = sqrt(g*g)
		g = (1.00/gn) * g
		push!(basis, g)
		next!(p)
	end

	serialize("basis$(dim)_$(M).jls", basis)
	return basis
end

function max_absval(f::MultiLinearFunction)
	return maximum([Base.abs(f.y[k]) for k in keys(f.y)])
end


function max_absval(f::PWMLinearFunction)
	return maximum([max_absval(f.pieces[k]) for k in keys(f.pieces)])
end


function constant_e1(D, gamma, basis)
	summand = 0.0
	for b in basis
		f_up = max_absval(b)
		summand += ((8 + 4*(2-gamma))/((gamma^2)))*(f_up^2)
	end
	println("summand = $summand")
	return summand * ((1/D) + (1/D^2)) 
end

function constant_sigma_2(D, gamma, basis)
	summand = 0.0
	for b in basis
		f_up = max_absval(b)
		summand += (((16*(f_up^4)))/(D^4))*((D * ((24- 36*gamma + 14*gamma^2 - gamma^3)/(gamma^4))) + ((6*D^2)*((4 - 4*gamma+ gamma^2)/(gamma^4))) + ((4*D^2)*((6 - 6*gamma + gamma^2)/(gamma^4))) )
	end
	summand *= length(basis)
	return summand
end

function cheby(D, basis, gamma, epsilon)
	e1 = constant_e1(D, gamma, basis)
	s2 = constant_sigma_2(D, gamma, basis)
	rad = e1 + epsilon
	delta = (s2/(epsilon^2))
	return rad, delta
end


function approximate_function(basis, f)
	approx = 0.0*basis[1] 
	coeffs = Vector{Float64}(undef, length(basis)) 
	
	# CAN MULTIPROCESS here using pmap to get coeffs
	if nworkers() > 1 && length(basis) > nworkers()
		s = ceil(Int, length(basis) / nworkers())
		chunks = [(s*m+1, basis[s*m + 1 : min(s*(m+1), length(basis))]) for m in 0:nworkers() - 1]
		result = @showprogress pmap(chunks) do basissubsettup
			start = basissubsettup[1]
			basissubset = basissubsettup[2]
			fs = 0.0*basis[1] 
			coeffs_inner = Dict{Int, Float64}()
			ind = 1
			for b in basissubset
				c = (b*f)[1]
				fs = fs + c*b	
				coeffs_inner[start + ind-1] = c
				ind += 1
			end
			return fs, coeffs_inner
		end
		merged = merge([r[2] for r in result]...)
		coeffs = [merged[i] for i in 1:length(basis)]
		for r in result
			approx = approx + r[1]
		end
	else
		ind = 1
		for b in basis
			println(ind/length(basis))
			c = (b*f)[1]
			approx = approx + c*b	
			coeffs[ind] = c 
			ind += 1
		end
	end
	return approx, coeffs
end

function approximate_kern(basis, data, stat)
	approx = 0.0*basis[1] 
	coeffs = Vector{Float64}(undef, length(basis)) 
	sta = x -> Func.eval_f(stat, [x]) 	
	# CAN MULTIPROCESS here using pmap to get coeffs
	if nworkers() > 1 && length(basis) > nworkers()
		s = ceil(Int, length(basis) / nworkers())
		chunks = [(s*m+1, basis[s*m + 1 : min(s*(m+1), length(basis))]) for m in 0:nworkers() - 1]
		result = @showprogress pmap(chunks) do basissubsettup
			start = basissubsettup[1]
			basissubset = basissubsettup[2]
			fs = 0.0*basis[1] 
			ind = 1
			tr = 1
			coeffs_inner = Dict{Int, Float64}()
			for b in basissubset
				println(tr / length(basissubset))
				tr += 1
				c = 0
				for d in data
					c += eval_f(b, d) / (sta(d[1]) *length(data))
				end
				coeffs_inner[start + tr-2] = c
				fs = fs + c*b	
				ind += 1
			end
			return fs, coeffs_inner
		end
		merged = merge([r[2] for r in result]...)
		coeffs = [merged[i] for i in 1:length(basis)]
		for r in result
			approx = approx + r[1]
		end
	else
		tr = 1
		for b in basis
			println(tr / length(basis))
			tr += 1
			c = 0.0
			for d in data
				c += eval_f(b, d) / (sta(d[1]) *length(data))
			end
			approx = approx + c*b	
			coeffs[tr - 1] = c
		end
	end
	return approx, coeffs
end

function approximate_emp(basis, data)
	approx = 0.0*basis[1] 
	coeffs = Vector{Float64}(undef, length(basis)) 
	if nworkers() > 1 && length(basis) > nworkers()
		s = ceil(Int, length(basis) / nworkers())
		chunks = [(s*m+1, basis[s*m + 1 : min(s*(m+1), length(basis))]) for m in 0:nworkers() - 1]
		result = @showprogress pmap(chunks) do basissubsettup
			start = basissubsettup[1]
			basissubset = basissubsettup[2]
			fs = 0.0*basis[1] 
			coeffs_inner = Dict{Int, Float64}()
			ind = 1
			rt = 1 
			for b in basissubset
				println(rt/length(basissubset))
				rt += 1
				c = 0.0
				for d in data
					c += Func.eval_f(b, d) / length(data)
				end
				coeffs_inner[start + rt-2] = c
				fs = fs + c*b	
				ind += 1
			end
			return fs, coeffs_inner
		end
		merged = merge([r[2] for r in result]...)
		coeffs = [merged[i] for i in 1:length(basis)]
		for r in result
			approx = approx + r[1]
		end
	else
		tr = 1
		D = length(data)
		for b in basis
			println(tr / length(basis))
			tr += 1
			c = 0
			@showprogress for d in data
				c += Func.eval_f(b, d) / D 
			end
			coeffs[tr - 1] = c
			approx = approx + c*b	
		end
	end
	
	return approx, coeffs
end


end # module

using .Func



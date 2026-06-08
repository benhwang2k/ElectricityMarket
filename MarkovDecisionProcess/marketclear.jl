module Market
using JuMP, Gurobi, CSV, DataFrames, Statistics
export PowerPlant, solve_market #, get_profit



# Might be interesting to also add minimum uptime, 
# but this gets more complicated.


# Maybe we should have devices modeled: IE controllable gen (ramps, capacity, price etc.)
# As well as renewable (peak capacity, stochasticity, cost etc).  not sure.
mutable struct PowerPlant
    name::String # fixed
    soc:Float64 # STATE VARIABLE units: kWH
    capacity::Float64 # fixed 
    min_output::Float64 # fixed
    no_load_cost::Float64 # fixed
    startup_cost::Float64 # fixed
    ramp_up::Float64 # fixed
    ramp_down::Float64 # fixed
    init_commit::Int
    init_gen::Float64
    block_fracs::Vector{Float64} # fixed
    block_costs::Vector{Float64} # fixed
    bid::Float64
end

function solve_market(plants, demands)
    model = Model(Gurobi.Optimizer)
    set_silent(model)

    I = eachindex(plants)
    H = eachindex(demands)

    #setting up variables g[i, h], u[i, h], y[i, h] (u, y are binary)
    @variable(model, generation[i in I, h in H] >= 0)
    @variable(model, commitment[i in I, h in H], Bin)
    @variable(model, startup[i in I, h in H], Bin)
    #setting up block generation
    @variable(model, block_generation[i in I, h in H, b in 1:length(plants[i].block_fracs)] >= 0)

    #minimize the cost, considering startup costs and strategic bidders
    @expression(model, startup_cost,
        sum(startup[i, h] * plants[i].startup_cost for i in I for h in H))
    @expression(model, no_load_cost, 
        sum(commitment[i, h] * plants[i].no_load_cost for i in I for h in H))
    @expression(model, block_variable_cost,
	    sum(block_generation[i, h, b] * plants[i].block_costs[b] * plants[i].bid 
        for i in I for h in H for b in 1:length(plants[i].block_fracs)))
    @objective(model, Min, startup_cost + no_load_cost + block_variable_cost)

    #upper bound, generation must be less than capacity
    @constraint(model, [i in I, h in H], generation[i, h] <= plants[i].capacity * commitment[i, h])

    #lower bound, if generator is on must produce min output
    @constraint(model, [i in I, h in H], generation[i, h] >= plants[i].min_output * commitment[i, h])

    #making it so startup cost is 1 when u goes from 0 to 1
    @constraint(model, [i in I, h in H[2:end]], startup[i, h] >= commitment[i, h] - commitment[i, h-1])

    #connecting block and normal gen, and capping each block
    @constraint(model, [i in I, h in H], generation[i, h] == 
        sum(block_generation[i, h, b] for b in 1:length(plants[i].block_fracs))
    )
    @constraint(model, [i in I, h in H, b in 1:length(plants[i].block_fracs)], 
        block_generation[i, h, b] <= plants[i].capacity * plants[i].block_fracs[b] * commitment[i, h]
    )


    #maybe include the startup edge cases, where ramp_up and ramp_down don't apply in that window, allowing for
    #generators to actaully start up if the min output is greater than their ramp up
    # Ben: hmmm. that might be nice, but then we also have to worry about minimum up-time calculations and it all gets a bit more hairy. IDK how much easier is it if we just assume that ramp up >> min output?
    #ramping up and down constraints on generation
    @constraint(model, [i in I, h in H[2:end]], generation[i, h] - generation[i, h-1] <= plants[i].ramp_up)
    @constraint(model, [i in I, h in H[2:end]], generation[i, h-1] - generation[i, h] <= plants[i].ramp_down)

    #initial condition constraints (so we can start from any combination of market horizons)
    @constraint(model, [i in I], generation[i, 1] - plants[i].init_gen <= plants[i].ramp_up)
    @constraint(model, [i in I], plants[i].init_gen - generation[i, 1] <= plants[i].ramp_down)
    @constraint(model, [i in I], startup[i, 1] >= commitment[i,1] -  plants[i].init_commit)


    #fufilling demand
    @constraint(model, demand_balance[h in H], sum(generation[i, h] for i in I) == demands[h])

    optimize!(model)

    undo = fix_discrete_variables(model)

    optimize!(model)

    price = dual.(demand_balance)
    gen_values = value.(generation)
    startup_values = value.(startup)
    commit_values = value.(commitment)
    block_gen = value.(block_generation)

    undo()

    profits = zeros(length(plants))

    for i in I
        for h in H
            variable_cost = sum(
                block_gen[i, h, b] * plants[i].block_costs[b]
                for b in 1:length(plants[i].block_fracs)
            )

            revenue = price[h] * gen_values[i, h]

            fixed_cost =
                startup_values[i, h] * plants[i].startup_cost +
                commit_values[i, h] * plants[i].no_load_cost

            profits[i] += revenue - variable_cost - fixed_cost
        end
    end

    return gen_values, price, profits
end

#im not sure but wouldnt it make more sense to have the profit function in the solve market?
#the solve market can just spit out the profits for each plant

#=
function get_profit(plant, price, block_gen, startup_values)
    profit = 0.0
    for h in H
	    #calculating profits from dual of demand constraint
	    profit += sum(price[h] - sum(block_gen[i, h, b] * plant.block_costs[b]))
	end return profit
end
=#
end

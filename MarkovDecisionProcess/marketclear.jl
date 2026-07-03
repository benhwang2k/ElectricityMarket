module Market
using JuMP, Gurobi, CSV, DataFrames, Statistics, MathOptInterface
export PowerPlant, solve_market #, get_profit

const MOI = MathOptInterface

mutable struct PowerPlant
    name::String
    capacity::Float64
    min_output::Float64
    no_load_cost::Float64
    startup_cost::Float64
    shutdown_cost::Float64
    min_up::Int
    min_down::Int
    ramp_up::Float64
    ramp_down::Float64
    init_commit::Int
    init_status::Int
    init_gen::Float64
    block_fracs::Vector{Float64}
    block_costs::Vector{Float64}
    bid::Float64
end

function solve_market(plants, demands)
    model = Model(Gurobi.Optimizer)
    set_silent(model)

    I = eachindex(plants)
    H = eachindex(demands)

    #setting up variables (u is commitment)
    @variable(model, generation[i in I, h in H] >= 0)
    @variable(model, commitment[i in I, h in H], Bin)
    @variable(model, startup[i in I, h in H], Bin)
    @variable(model, shutdown[i in I, h in H], Bin)

    #setting up block generation
    @variable(model, block_generation[i in I, h in H, b in 1:length(plants[i].block_fracs)] >= 0)

    #minimize the cost, considering startup costs and strategic bidders
    @expression(model, startup_cost,
        sum(startup[i, h] * plants[i].startup_cost for i in I for h in H))
    @expression(model, shutdown_cost,
        sum(shutdown[i, h] * plants[i].shutdown_cost for i in I for h in H))
    @expression(model, no_load_cost, 
        sum(commitment[i, h] * plants[i].no_load_cost for i in I for h in H))
    @expression(model, block_variable_cost,
	    sum(block_generation[i, h, b] * plants[i].block_costs[b] * plants[i].bid 
        for i in I for h in H for b in 1:length(plants[i].block_fracs)))
    @objective(model, Min, startup_cost + shutdown_cost + no_load_cost + block_variable_cost)

    #upper bound, generation must be less than capacity
    @constraint(model, [i in I, h in H], generation[i, h] <= plants[i].capacity * commitment[i, h])

    #lower bound, if generator is on must produce min output
    @constraint(model, [i in I, h in H], generation[i, h] >= plants[i].min_output * commitment[i, h])

    #creating startup and shutdown binary variables
    @constraint(model, [i in I], commitment[i, 1] - plants[i].init_commit == startup[i, 1] - shutdown[i, 1])
    @constraint(model, [i in I, h in H[2:end]], commitment[i, h] - commitment[i, h-1] == startup[i, h] - shutdown[i, h])

    #connecting block and normal gen, and capping each block
    @constraint(model, [i in I, h in H], generation[i, h] == 
        sum(block_generation[i, h, b] for b in 1:length(plants[i].block_fracs))
    )
    @constraint(model, [i in I, h in H, b in 1:length(plants[i].block_fracs)], 
        block_generation[i, h, b] <= plants[i].capacity * plants[i].block_fracs[b] * commitment[i, h]
    )

    #ramping up and down constraints on generation
    @constraint(model, [i in I, h in H[2:end]], generation[i, h] - generation[i, h-1] <= plants[i].ramp_up)
    @constraint(model, [i in I, h in H[2:end]], generation[i, h-1] - generation[i, h] <= plants[i].ramp_down)

    #initial condition constraints (so we can start from any combination of market horizons)
    @constraint(model, [i in I], generation[i, 1] - plants[i].init_gen <= plants[i].ramp_up)
    @constraint(model, [i in I], plants[i].init_gen - generation[i, 1] <= plants[i].ramp_down)
    

    #minimum up and downtime constraints
    for i in I
        T = length(H)
        UT = plants[i].min_up
        G = min(T, max(0, (UT - plants[i].init_status) * plants[i].init_commit))
        # (21)
        @constraint(model, sum(1 - commitment[i, h] for h in 1:G) == 0)

        # (22)
        if max(G + 1, 2) <= T - UT + 1
            @constraint(model, [h in max(G + 1, 2):(T-UT+1)],
                sum(commitment[i, n] for n in h:(h+UT-1)) >= UT * (startup[i, h])
            )
        end

        # (23)
        if T - UT + 2 <= T
            @constraint(model, [h in (T-UT+2):T],
                sum(commitment[i, n] - (startup[i, h]) for n in h:T) >= 0
            )
        end
    end

    for i in I
        T = length(H)
        DT = plants[i].min_down
        L = min(T, max(0, (DT + plants[i].init_status) * (1 - plants[i].init_commit)))
        # (24)
        @constraint(model, sum(commitment[i, h] for h in 1:L) == 0)

        # (25)
        if max(L + 1, 2) <= T - DT + 1
            @constraint(model, [h in max(L + 1, 2):(T-DT+1)],
                sum(1 - commitment[i, n] for n in h:(h+DT-1)) >= DT * (shutdown[i, h])
            )
        end

        # (26)
        if T - DT + 2 <= T
            @constraint(model, [h in (T-DT+2):T],
                sum(1 - commitment[i, n] - (shutdown[i, h]) for n in h:T) >= 0
            )
        end
    end

    #fufilling demand
    @constraint(model, demand_balance[h in H], sum(generation[i, h] for i in I) == demands[h])

    optimize!(model)
    
    if termination_status(model) != MOI.OPTIMAL
        println("Status: ", termination_status(model))
        return nothing, nothing, nothing
    end

    undo = fix_discrete_variables(model)

    optimize!(model)

    price = dual.(demand_balance)
    gen_values = value.(generation)
    #startup_values = value.(startup)
    commit_values = value.(commitment)
    #block_gen = value.(block_generation)

    undo()
    
    #= old profits code
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
    =#

    return gen_values, price, commit_values
end

#im not sure but wouldnt it make more sense to have the profit function in the solve market?
#the solve market can just spit out the profits for each plant

#=
function get_profit(plant, price, block_gen, startup_values)
    profit = 0.0
    for h in H
	    #calculating profits from dual of demand constraint
	    profit += sum(price[h] - sum(block_gen[i, h, b] * plant.block_costs[b]))
	end
    return profit
end
=#
end

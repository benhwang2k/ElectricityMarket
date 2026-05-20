using JuMP, Gurobi, CSV, DataFrames, Statistics

struct PowerPlant
    name::String
    capacity::Float64
    min_output::Float64
    variable_cost::Float64
    startup_cost::Float64
    ramp_up::Float64
    ramp_down::Float64
    init_commit::Int
    init_gen::Float64
    is_strategic::Bool
end

struct Consumer 
    name::String
    capacity::Float64
    min_output::Vector{Float64}
    variable_cost::Vector{Float64}
end

plants_df = CSV.read("MarketClearing/data/plants.csv", DataFrame)
plants = PowerPlant[]
for i in eachindex(plants_df.name)
    push!(plants, PowerPlant(
        plants_df.name[i],
        plants_df.capacity[i],
        plants_df.min_output[i],
        plants_df.variable_cost[i],
        plants_df.startup_cost[i],
        plants_df.ramp_up[i],
        plants_df.ramp_down[i],
        plants_df.init_commit[i],
        plants_df.init_gen[i],
        Bool(plants_df.is_strategic[i])
    ))
end

consumers_df = CSV.read("MarketClearing/data/consumers.csv", DataFrame)
consumers = Consumer[]
for i in 1:size(consumers_df)[1] 
    push!(consumers, Consumer(consumer_df[1,1], consumer_df[1,2], consumer_df[1,3:(3+T)], consumer_df[1,(3+T+1):(3+2*T+1)]))
end



function solve_market(plants, demands, k)
    model = Model(Gurobi.Optimizer)
    set_silent(model)

    I = eachindex(plants)
    H = eachindex(demands)

    #setting up variables g[i, h], u[i, h], y[i, h] (u, y are binary)
    @variable(model, generation[i in I, h in H] >= 0)
    @variable(model, commitment[i in I, h in H], Bin)
    @variable(model, startup[i in I, h in H], Bin)

    #minimize the cost, considering startup costs and strategic bidders
    @expression(model, startup_cost,
        sum(startup[i, h] * plants[i].startup_cost for i in I for h in H))
    @expression(model, variable_cost,
        sum(generation[i, h] * (plants[i].is_strategic ? k * plants[i].variable_cost : 
        plants[i].variable_cost) for i in I for h in H))
    @objective(model, Min, startup_cost + variable_cost)

    #upper bound, generation must be less than capacity
    @constraint(model, [i in I, h in H], generation[i, h] <= plants[i].capacity * commitment[i, h])

    #lower bound, if generator is on must produce min output
    @constraint(model, [i in I, h in H], generation[i, h] >= plants[i].min_output * commitment[i, h])


    #making it so startup cost is 1 when u goes from 0 to 1
    @constraint(model, [i in I, h in H[2:end]], startup[i, h] >= commitment[i, h] - commitment[i, h-1])

    #maybe include the startup edge cases, where ramp_up and ramp_down don't apply in that window, allowing for
    #generators to actaully start up if the min output is greater than their ramp up

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
    
    gen_values = [value(generation[i, h]) for i in I, h in H]
    startup_values = [value(startup[i,h]) for i in I, h in H]
    
    #fixing commitment values and then rerunning the program
    undo = fix_discrete_variables(model);

    optimize!(model)
    price = dual.(demand_balance)

    undo()

    profit = 0.0
    strategic_generation = 0.0

    for i in I
        for h in H
            if plants[i].is_strategic
                strategic_generation += gen_values[i, h]
                #calculating profits from dual of demand constraint
                profit += (price[h] - plants[i].variable_cost) * gen_values[i, h] - 
                    startup_values[i, h] * plants[i].startup_cost
            end
        end
    end
    return gen_values, price, profit, strategic_generation
end

function get_profit(plant, gen_values, startup_values)
    profit = 0.0
    strategic_generation = 0.0

    for h in H
	if plant.is_strategic
	    strategic_generation += gen_values[i, h]
	    #calculating profits from dual of demand constraint
	    profit += (price[h] - plants[i].variable_cost) * gen_values[i, h] - 
		startup_values[i, h] * plants[i].startup_cost
	end
    end
    return profit
end

results = DataFrame(
    k = Float64[],
    avg_price = Float64[],
    profit = Float64[],
    strategic_generation = Float64[]
)


for k in 1.0:0.1:2.0
    dispatch, price, profit, strategic_gen = solve_market(plants, demands_df.demand, k)

    push!(results, (
        k = k,
        avg_price = mean(price),
        profit = profit,
        strategic_generation = strategic_gen
    ))
end 
println(results)
CSV.write("MarketClearing/data/results.csv", results)

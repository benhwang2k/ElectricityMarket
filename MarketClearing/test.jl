using JuMP, Gurobi, CSV, DataFrames, Statistics

struct PowerPlant
    name::String
    capacity::Float64
    min_output::Float64
    variable_cost::Float64
    startup_cost::Float64
    is_strategic::Bool
end

plants_df = CSV.read("03_data/plants.csv", DataFrame)
plants = PowerPlant[]
for i in eachindex(plants_df.name)
    push!(plants, PowerPlant(
        plants_df.name[i],
        plants_df.capacity[i],
        plants_df.min_output[i],
        plants_df.variable_cost[i],
        plants_df.startup_cost[i],
        Bool(plants_df.is_strategic[i])
    ))
end

demands_df = CSV.read("03_data/demand.csv", DataFrame)
demands = demands_df.demand

model = Model(Gurobi.Optimizer)
#set_silent(model)

I = eachindex(plants)
H = eachindex(demands)

#setting up variables g[i, h] and u[i, h] with u being binary
@variable(model, generation[i in I, h in H] >= 0)
@variable(model, commitment[i in I, h in H], Bin)

#minimize the cost, considering startup costs and strategic bidders
@expression(model, startup_cost,
    sum(commitment[i, h] * plants[i].startup_cost for i in I for h in H))
@expression(model, variable_cost,
    sum(generation[i, h] * (plants[i].is_strategic ? 1.5 * plants[i].variable_cost : 
    plants[i].variable_cost) for i in I for h in H))
@objective(model, Min, startup_cost + variable_cost)

#upper bound, generation must be less than capacity
@constraint(model, [i in I, h in H], generation[i, h] <= plants[i].capacity * commitment[i, h])

#lower bound, if generator is on must produce min output
@constraint(model, [i in I, h in H], generation[i, h] >= plants[i].min_output * commitment[i, h])

#fufilling demand
@constraint(model, demand_balance[h in H], sum(generation[i, h] for i in I) == demands[h])

grb = backend(model)

optimize!(model)
    
MOI.get(grb, Gurobi.ModelAttribute("Status"))

println()
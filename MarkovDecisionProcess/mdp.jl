using .Market
using CSV, DataFrames
using Random
###
# Add state, action, environmnet, reward 
#
# state: locally cleared generation and price, also capacity and other internal variables
# action: bid
# environment: market clearing <- lower level opt, and competitor's bids
# reward: profit 
# 		eventually, this can be calculated by a separate
# 		tracking problem, but for now lets just assume we have
# 		perfect tracking capability


# Each real time market is calculated for a number of periods T
T = 1


# Cleared generation for a horizon T
g = [0.0]*T
# Price at each time
l = [0.0]*T
# commitment for each time
u = [false]*T


mutable struct State
  gen::Matrix{Float64}
  price::Vector{Float64}
end


struct Bid
  competition::Vector{Float64}
end


plants_df = CSV.read("MarkovDecisionProcess/data/plants.csv", DataFrame)
plants = Market.PowerPlant[]
for i in eachindex(plants_df.name)
    push!(plants, Market.PowerPlant(
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


demands_df = CSV.read("MarketClearing/data/demand.csv", DataFrame)


function next_state(plants, state, new_demand)
  if !(is_demand_feasible(plants, new_demand))
    print("demand not feasible")
    return state
  end
  # capactiy of next day = cap of prev - gen + solar
  k = 1
  
  # solve
  g, price, profit, one_strategic_gen = solve_market(plants, new_demand, k)
  next_state = State(g, price)  

  # plant controllable capacity (ie gas generation)
  # plus a stochastic renewable resource
  # storage term
  gas_gen = [200, 300, 100, 250]
  renew_gen = [100, 100, 100, 100]
  battery_cap = [100, 100, 100, 100]
  for i in 1:4 
    plants[i].capacity = max(min(battery_cap[i], plants[i].capacity - gas_gen[i] -state.gen[i]),0) + renew_gen[i]*rand() + gas_gen[i] 
    if state.gen[i] > 0
      plants[i].init_commit = 1 
    else
      plants[i].init_commit = 0
    end
    plants[i].init_gen = state.gen[i] 
  end
  return next_state
end


function is_demand_feasible(plants, demand)
  downcapacity = 0
  capacity = 0
  for p in plants
    upc = min(p.init_gen + p.ramp_up, p.capacity)
    downc = max(p.init_gen - p.ramp_down, p.min_output)
    if p.init_gen < p.ramp_down
      downc = 0
    end
    capacity += upc
    downcapacity += downc
  end
  return downcapacity <= demand <= capacity
end


function print_cap(plants)
  println()
  println()
  println("capacities") 
  println([p.capacity for p in plants])
  println()
  println()
end


function simulate()
  plants = Market.PowerPlant[]
  for i in eachindex(plants_df.name)
      push!(plants, Market.PowerPlant(
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

  # constant load
  demand = 300 
  # solve
  g, price, profit, one_strategic_gen = solve_market(plants, demand, 1)
  s = State(g, price)  
  state_array = [s]
  for t in 1:5
    if t > 2
      demand = 400
    end
    print("t = ")
    print(t)
    println()
    s = next_state(plants, s, demand)
    push!(state_array, s)
    print_cap(plants)
  end
  return state_array, plants
end






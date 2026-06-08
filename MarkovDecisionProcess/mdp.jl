using Pkg 
Pkg.activate(joinpath(@__DIR__, ".."))

include("marketclear.jl")
using .Market
using CSV, DataFrames
using Random
###
# Add state, action, environment, reward 
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
  plants::Vector{PowerPlant}
end


struct Bid
  competition::Vector{Float64}
end


plants_df = CSV.read("MarkovDecisionProcess/data/plants.csv", DataFrame)
function load_plants()
  plants = Market.PowerPlant[]
  for i in eachindex(plants_df.name)
    block_fracs = [
      plants_df.block1_frac[i],
      plants_df.block2_frac[i],
      plants_df.block3_frac[i],
    ]

    block_costs = [
      plants_df.block1_cost[i],
      plants_df.block2_cost[i],
      plants_df.block3_cost[i],
    ]

    if !isapprox(sum(block_fracs), 1.0; atol=1e-6)
      error("Block fractions for $(plants_df.name[i]) must sum to 1.")
    end

    push!(plants, Market.PowerPlant(
      plants_df.name[i],
      0.0, #state of charge starts at zero
      plants_df.capacity[i],
      plants_df.min_output[i],
      plants_df.no_load_cost[i],
      plants_df.startup_cost[i],
      plants_df.ramp_up[i],
      plants_df.ramp_down[i],
      plants_df.init_commit[i],
      plants_df.init_gen[i],
      block_fracs,
      block_costs,
      plants_df.bid[i])
    )
  end
  return plants
end


demands_df = CSV.read("MarketClearing/data/demand.csv", DataFrame)


function dispatch(plant, gen)
  """
  allocate the controllable, renewable, and battery resources
  to provide generation = gen. Return the cost of doing so.

  Also update the plant state variables
  """
  
  return 1 
end

function real_time_execution(plant_ind, plants, g, price)
  """
  This function uses the plant to track the 
  generation profile (g). The result is the 
  profile of the plant's capacity as a function
  of time as well as the profit of the plant.
  """
  
  # Control law says how to choose between battery and generator
  p = plants[plant_ind]
  T = size(g)[2]

  reward = 0 
  for t in 1:T
    # plant will supply the generation.-- update plant state and return reward
    reward += dispatch(p, g[plant_ind, t])
  end
  return reward
end



# This function models both the market cleared, and the real-time execution 
# of the cleared generation.
function next_state(state, new_demand)
  plants = state.plants
  if !(is_demand_feasible(plants, new_demand))
    print("demand not feasible")
    return state
  end
  
  # Upper Level problem: optimize bid for day ahead market.
  k = 1
  
  # Lower Level Problem: Market Clearing! 
  g, price = Market.solve_market(plants, [new_demand])

  # Real-Time Execution
  # Calculate the plant's real-time market capacity
  rewards = []
  for i in 1:4 
    push!(rewards, sum(price[:].* g[i,:]) - real_time_execution(i, plants, g, price))
  end
  # Record the previous market data as a 'state' to transition from.
  next_state = State(g, price, plants)  

  # update the commitment and initial values of the plants for the next market 
  for i in 1:4 
    if state.gen[i] > 0
      plants[i].init_commit = 1 
    else
      plants[i].init_commit = 0
    end
    plants[i].init_gen = state.gen[i] 
  end
  return next_state, rewards
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
  plants = load_plants()

  # constant load
  demand = 300 
  # solve
  g, price, profits = Market.solve_market(plants, [demand])
  s = State(g, price)  
  state_array = [s]
  for t in 1:5
    if t > 2
      demand = 400
    end
    print("t = ")
    print(t)
    println()
    s = next_state(s, demand)
    push!(state_array, s)
    print_cap(plants)
  end
  return state_array, plants
end






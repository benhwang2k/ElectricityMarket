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
end


struct Bid
  competition::Vector{Float64}
end


plants_df = CSV.read("MarkovDecisionProcess/data/plants.csv", DataFrame)
function load_plants()
  plants = Market.PowerPlant[]
  for i in eachindex(plants_df.name)
    init_commit = 0
    if (plants_df.init_status[i] > 0)
      init_commit = 1
    end

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
      plants_df.capacity[i],
      plants_df.min_output[i],
      plants_df.no_load_cost[i],
      plants_df.startup_cost[i],
      plants_df.shutdown_cost[i],
      plants_df.min_up[i],
      plants_df.min_down[i],
      plants_df.ramp_up[i],
      plants_df.ramp_down[i],
      init_commit,
      plants_df.init_status[i],
      plants_df.init_gen[i],
      block_fracs,
      block_costs,
      plants_df.bid[i])
    )
  end
  return plants
end


demands_df = CSV.read("MarketClearing/data/demand.csv", DataFrame)


function next_state(plants, new_demand)
  if !(is_demand_feasible(plants, new_demand))
    print("demand not feasible")
    return state
  end
  # capactiy of next day = cap of prev - gen + solar
  k = 1
  
  # solve
  g, price, commitment = Market.solve_market(plants, [new_demand])
  

  if g === nothing
    println("market solve failed")
    return nothing
  end

  next_state = State(g, price)  
  # plant controllable capacity (ie gas generation)
  # plus a stochastic renewable resource
  # storage term
  gas_gen = [300, 200, 100, 250]
  renew_gen = [100, 100, 100, 100]
  battery_cap = [100, 100, 100, 100]
  for i in eachindex(plants)
    current_gen = g[i, 1]
    current_commit = round(Int, commitment[i, 1])

    old_commit = plants[i].init_commit
    old_status = plants[i].init_status

    #update capacity using current gen
    plants[i].capacity = max(min(battery_cap[i], plants[i].capacity - gas_gen[i] - current_gen), 0.0) + renew_gen[i] * rand() + gas_gen[i]

    #update commitment and gen initial conditioon
    plants[i].init_commit = current_commit
    plants[i].init_gen = current_gen

    if current_commit == 1
      if old_commit == 1
        plants[i].init_status = max(1, old_status + 1)
      else
        plants[i].init_status = 1
      end
    else
      if old_commit == 0
        plants[i].init_status = min(-1, old_status - 1)
      else
        plants[i].init_status = -1
      end
    end
  end
  #=
  for i in 1:4 
    plants[i].capacity = max(min(battery_cap[i], plants[i].capacity - gas_gen[i] -state.gen[i]),0) + renew_gen[i]*rand() + gas_gen[i] 
    if state.gen[i] > 0
      plants[i].init_commit = 1 
    else
      plants[i].init_commit = 0
    end
    plants[i].init_gen = state.gen[i] 
  end =#
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
  plants = load_plants()
  state_array = State[]

  for t in 1:5
    demand = t > 2 ? 400.0 : 300.0

    println("t = ", t)

    s = next_state(plants, demand)

    if s === nothing
        break
    end

    push!(state_array, s)
    print_cap(plants)
  end

  return state_array, plants
end






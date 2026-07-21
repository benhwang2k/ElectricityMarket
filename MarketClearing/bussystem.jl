using JuMP
using Gurobi
import MathOptInterface as MOI

struct Bus
    id::Int
    p_demand::Float64  #active demand, MW
    q_demand::Float64  #reactive demand, MVAr
end

struct Branch
    from_bus::Int
    to_bus::Int
    r::Float64     #resistance, per unit
    x::Float64     #reactance, per unit
    rate::Float64  #in MVA
end

struct PowerPlant
    name::String
    bus::Int
    p_max::Float64
    q_min::Float64
    q_max::Float64
    cost::Float64
end

buses = [
    Bus(1, 0.0, 0.0),
    Bus(2, 100.0, 30.0),
    Bus(3, 50.0, 15.0),
]

branches = [
    Branch(1, 2, 0.01, 0.10, 100.0),
    Branch(1, 3, 0.02, 0.20, 100.0),
    Branch(2, 3, 0.025, 0.25, 60.0),
]

plants = [
    PowerPlant("Cheap Plant", 1, 200.0, -100.0, 100.0, 20.0),
    PowerPlant("Expensive Plant", 3, 200.0, -100.0, 100.0, 60.0),
]

function solve_bus_market(plants, buses, branches; baseMVA = 100.0, v_min = 0.95, v_max = 1.05)
    model = Model(Gurobi.Optimizer)
    set_silent(model)
    set_optimizer_attribute(model, "QCPDual", 1)

    I = eachindex(plants)
    V = [bus.id for bus in buses]
    E = eachindex(branches)

    @assert length(unique(V)) == length(V) "Bus IDs must be unique."
    @assert all(
        branch.from_bus in V && branch.to_bus in V
        for branch in branches
    ) "Every branch endpoint must correspond to a bus."
    @assert all(
        branch.r >= 0 &&
        branch.x > 0 &&
        branch.rate > 0
        for branch in branches
    ) "Branch parameters must be valid."

    #converting to per unit
    p_demand_pu = Dict(
        bus.id => bus.p_demand / baseMVA
        for bus in buses
    )

    q_demand_pu = Dict(
        bus.id => bus.q_demand / baseMVA
        for bus in buses
    )

    rate_pu = Dict(
        e => branches[e].rate / baseMVA
        for e in E
    )

    # Variables
    @variable(model, 0 <= p_generation[i in I] <= plants[i].p_max / baseMVA)
    @variable(model, plants[i].q_min / baseMVA <= q_generation[i in I] <= plants[i].q_max / baseMVA)

    #voltage variables are squared so limits are squared
    @variable(model, v_min^2 <= voltage_sq[v in V] <= v_max^2)
    @variable(model, current_sq[e in E] >= 0)
    @variable(model, active_flow[e in E])
    @variable(model, reactive_flow[e in E])

    # Fix voltage magnitude at the reference bus
    reference_bus = V[1]
    @constraint(model, voltage_sq[reference_bus] == 1.0)

    #voltage change across branches (1a)
    @constraint(model, [e in E], voltage_sq[branches[e].to_bus] == voltage_sq[branches[e].from_bus] -
        (2 * (branches[e].r * active_flow[e] + branches[e].x * reactive_flow[e])) + 
        (branches[e].r^2 + branches[e].x^2) * current_sq[e]
    )

    #Power current relationship - relaxed with inequality (1d)
    @constraint(model, [e in E], active_flow[e]^2 + reactive_flow[e]^2 <= 
        voltage_sq[branches[e].from_bus] * current_sq[e]
    )

    #sending end apparent power limit
    @constraint(model, [e in E],
        active_flow[e]^2 + reactive_flow[e]^2 <=
        rate_pu[e]^2
    )

    #receiving end apparent power limit
    @constraint(model, [e in E],
        (active_flow[e] - branches[e].r * current_sq[e])^2 +
        (reactive_flow[e] - branches[e].x * current_sq[e])^2 <=
        (rate_pu[e])^2
    )
    #------------------------------------------------------------------------
    #branch flow model 
    p_balance_expr = Dict(
        v => AffExpr(0.0)
        for v in V
    )

    q_balance_expr = Dict(
        v => AffExpr(0.0)
        for v in V
    )

    for i in I
        plant_bus = plants[i].bus

        add_to_expression!(
            p_balance_expr[plant_bus],
            1.0,
            p_generation[i]
        )

        add_to_expression!(
            q_balance_expr[plant_bus],
            1.0,
            q_generation[i]
        )
    end


    for e in E
        branch = branches[e]

        from_bus = branch.from_bus
        to_bus = branch.to_bus

        # Power leaves the from-bus
        add_to_expression!(
            p_balance_expr[from_bus],
            -1.0,
            active_flow[e]
        )

        add_to_expression!(
            q_balance_expr[from_bus],
            -1.0,
            reactive_flow[e]
        )

        # Active power arrives at the to-bus:
        # P_ij - r_ij * J_ij
        add_to_expression!(
            p_balance_expr[to_bus],
            1.0,
            active_flow[e]
        )

        add_to_expression!(
            p_balance_expr[to_bus],
            -branch.r,
            current_sq[e]
        )

        # Reactive power arrives at the to-bus:
        # Q_ij - x_ij * J_ij
        add_to_expression!(
            q_balance_expr[to_bus],
            1.0,
            reactive_flow[e]
        )

        add_to_expression!(
            q_balance_expr[to_bus],
            -branch.x,
            current_sq[e]
        )
    end

    active_balance = Dict{Int, ConstraintRef}()
    reactive_balance = Dict{Int, ConstraintRef}()

    #equating demands to the power
    for v in V
        active_balance[v] = @constraint(
            model,
            p_balance_expr[v] == p_demand_pu[v]
        )

        reactive_balance[v] = @constraint(
            model,
            q_balance_expr[v] == q_demand_pu[v]
        )
    end
    #------------------------------------------------------------------------
    #minimize the cost
    @objective(model, Min, 
        sum(plants[i].cost * baseMVA * p_generation[i] for i in I)
    )

    optimize!(model)

    status = termination_status(model)

    if status != MOI.OPTIMAL
        error("Model did not solve to optimality. Status: $status")
    end

    #collecting our values
    p_generation_values = Dict(i => baseMVA * value(p_generation[i]) for i in I)
    q_generation_values = Dict(i => baseMVA * value(q_generation[i]) for i in I)
    voltage_values = Dict(v => sqrt(value(voltage_sq[v])) for v in V)
    current_sq_values = Dict(e => value(current_sq[e]) for e in E)
    active_flow_values = Dict(e => baseMVA * value(active_flow[e]) for e in E)
    reactive_flow_values = Dict(e => baseMVA * value(reactive_flow[e]) for e in E)

    #computing important info from values
    branch_values = Dict{Int, Any}()

    for e in E
        branch = branches[e]

        p_send_pu = value(active_flow[e])
        q_send_pu = value(reactive_flow[e])
        j_value = value(current_sq[e])

        p_receive_pu =
            p_send_pu - branch.r * j_value

        q_receive_pu =
            q_send_pu - branch.x * j_value

        sending_mva =
            baseMVA * hypot(p_send_pu, q_send_pu)

        receiving_mva =
            baseMVA * hypot(p_receive_pu, q_receive_pu)

        branch_values[e] = (
            from_bus = branch.from_bus,
            to_bus = branch.to_bus,

            p_send_mw =
                baseMVA * p_send_pu,

            q_send_mvar =
                baseMVA * q_send_pu,

            p_receive_mw =
                baseMVA * p_receive_pu,

            q_receive_mvar =
                baseMVA * q_receive_pu,

            active_loss_mw =
                baseMVA * branch.r * j_value,

            reactive_loss_mvar =
                baseMVA * branch.x * j_value,

            sending_mva = sending_mva,
            receiving_mva = receiving_mva,

            rate_mva = branch.rate,

            sending_loading_percent =
                100 * sending_mva / branch.rate,

            receiving_loading_percent =
                100 * receiving_mva / branch.rate
        )
    end

    price_values = Dict(v => -shadow_price(active_balance[v]) / baseMVA for v in V)
    return (
        p_generation_mw = p_generation_values,
        q_generation_mvar = q_generation_values,
        voltage_pu = voltage_values,
        current_sq_pu = current_sq_values,
        branches = branch_values,
        prices = price_values,
        objective = objective_value(model)
    )

end

result = solve_bus_market(plants, buses, branches)

println("\nActive generation, MW:")
println(result.p_generation_mw)

println("\nReactive generation, MVAr:")
println(result.q_generation_mvar)

println("\nVoltage magnitude, p.u.:")
println(result.voltage_pu)

println("\nSquared current, p.u.:")
println(result.current_sq_pu)

println("\nBranches:")
for (e, values) in result.branches
    println("Branch $e: $values")
end

println("\nActive-power prices, \$/MWh:")
println(result.prices)

println("\nObjective:")
println(result.objective)
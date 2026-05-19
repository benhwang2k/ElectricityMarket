using JuMP, Gurobi


function build_model()
  model = Model(Gurobi.Optimizer)

  @variable(model, x >= 0)
  @variable(model, 0 <= y <= 3)
  @objective(model, Min, x + y)
  @constraint(model, c1, x + y >= 4)
  model
end

greeting() = "Hello earthling"

function sol()
  model = build_model()
  optimize!(model)
end


###
# My REPL uses the Revise package to update any
# code in the REPL automatically as along as the code 
# is inside a function
#
# This is how to format a Gurobi optimzer so that it makes
# development in the REPL easy without having to inlude the file 
# again and again
#
# using Revise
# includet("simple.jl")
# sol()
#
# ----- update simple.jl and save -----
#
# sol() # provides the updated output.
#

using Pkg
Pkg.activate(".")
Pkg.instantiate()

using MPSGE
using DataFrames
using PlotlyJS
import JuMP.Containers: DenseAxisArray, @container

import JuMP
using Ipopt


include("indexed_models.jl")

## Ensure calibration

d4_data = IndexedData();
d4_model = dynamic_index_model(d4_data);
solve!(d4_model, cumulative_iteration_limit=0)


## Set up forecasted quantities

goods = d4_data.goods
time_periods = d4_data.time_periods
X0 = d4_data.X


# We are going to assume a constant 2% growth rate for good A
# For good B, we will assume 2.2% growth until 2050 and then 1.6% thereafter

@container(SECQREF[g∈goods, t∈time_periods], 
    g == :A ? 1.02^(t - 2000) : (t<2051 ? 1.022^(t-2000) : 1.022^50 * 1.016^(t - 2050))
)

@container(AVGREF_corrected[t∈time_periods], 
    sum(SECQREF[g,t]*X0[g] for g in goods)/sum(X0[g] for g in goods)
)


## Fix the new QREF values, and use this to determine the values of SK

set_lower_bound.(d4_model[:SK], -.99)
set_value!.(d4_model[:QREFJ], SECQREF)
solve!(d4_model)

## Adjust the PHI values based on the new SK values. PHI = 1 + SK

set_value!.(d4_model[:PHI], 1 .+value.(d4_model[:SK]))
fix.(d4_model[:SK], 0)
fix.(d4_model[:ADJ_PKJ], 0)

# Check that we are out of equilibrium
solve!(d4_model, cumulative_iteration_limit=0)

benchmark_values = DataFrame(
    t = d4_data.time_periods,
    Y = [value(d4_model[:Y][t]) for t in d4_data.time_periods],
    C = [value(d4_model[:C][t]) for t in d4_data.time_periods],
    I = [value(d4_model[:I][t]) for t in d4_data.time_periods],
    K = [value(d4_model[:K][t]) for t in d4_data.time_periods],
    X1 = [value(d4_model[:X][:A, t]) for t in d4_data.time_periods],
    X2 = [value(d4_model[:X][:B, t]) for t in d4_data.time_periods],
) |>
x -> stack(x, Not(:t), variable_name = :variable, value_name = :benchmark)

# Set a shock and solve
set_value!.(d4_model[:TAX][2010:end], 0.1);
solve!(d4_model)

shock_values = DataFrame(
    t = d4_data.time_periods,
    Y = [value(d4_model[:Y][t]) for t in d4_data.time_periods],
    C = [value(d4_model[:C][t]) for t in d4_data.time_periods],
    I = [value(d4_model[:I][t]) for t in d4_data.time_periods],
    K = [value(d4_model[:K][t]) for t in d4_data.time_periods],
    X1 = [value(d4_model[:X][:A, t]) for t in d4_data.time_periods],
    X2 = [value(d4_model[:X][:B, t]) for t in d4_data.time_periods],
) |>
x -> stack(x, Not(:t), variable_name = :variable, value_name = :shock)



df = leftjoin(
    benchmark_values,
    shock_values,
    on = [:t, :variable]
) |>
x -> transform(x,
    [:benchmark, :shock] => ByRow((b,s) -> 100*(s/b - 1)) => :value
) |>
x -> stack(x, Not(:t, :variable), variable_name = :model, value_name = :value) 



## Graphs


mode_order = ["Y", "C", "I", "K", "X1", "X2"]
modes = Dict(
     "Y" => "Output",
     "C" => "Consumption",
     "I" => "Investment",
     "K" => "Capital",
     "X1" => "Output A",
     "X2" => "Output B",
)


layout = Layout(
    updatemenus = [
        attr(
            active=0,
            buttons = [
                attr(
                    label=modes[mode],
                    method="update",
                    args=[
                        attr(visible = [m == mode for m in mode_order] ),
                        attr(title=modes[mode])
                    ]
                ) for mode in mode_order
            ]
        )
    ];
    title = string(modes[mode_order[1]]),
    yaxis_title="% deviation from baseline",
    xaxis_title="Year"
)

df |>
    x -> plot(
        vec([
            scatter(
                df |> x-> subset(x, :variable => ByRow(==(mode)), :model => ByRow(==(m))),
                x=:t,
                y=:value,
                #line_color=:model,
                mode="lines",
                text=:model,
                name = string(m),
                visible = mode == mode_order[1]
            ) for mode in mode_order, m in unique(df.model)
        ]),
        layout
    )


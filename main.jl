using Pkg
Pkg.activate(".")
Pkg.instantiate()

using MPSGE
using DataFrames
using PlotlyJS
import JuMP.Containers: DenseAxisArray, @container


## Include model definitions
include("scalar_data.jl")
include("scalar_models.jl")





## Static Model ##

static_data = StaticData(
    200, # Y0
    180, # C0
    150, # LS0
    50,  # KS0
    20  # I0
    )

static_model = static_scalar_model(static_data)
solve!(static_model, cumulative_iteration_limit=0)
generate_report(static_model)

set_value!(static_model[:TAX], 0.1);
solve!(static_model)
generate_report(static_model)

#######################################
## Interest Rate Determined by Model ##
#######################################

d1_data = dynamic_data_1(static_data);
d1_model = dynamic_scalar_model(d1_data);
solve!(d1_model, cumulative_iteration_limit=0)
set_value!.(d1_model[:TAX][2010:2100], 0.1);
solve!(d1_model)
df1 = dyn_model_report(d1_model, d1_data; value_name = :Fixed_Rate)


############################################
## Interest Rate Given, adjust investment ##
############################################

d2_data = dynamic_data_2(static_data);
d2_model = dynamic_scalar_model(d2_data);
solve!(d2_model, cumulative_iteration_limit=0)
set_value!.(d2_model[:TAX][2010:2100], 0.1);
solve!(d2_model)
df2 = dyn_model_report(d2_model, d2_data; value_name = :Adjust_Investment)


############################################
## Interest Rate Given, adjust capital #####
############################################
d3_data = dynamic_data_3(static_data);
d3_model = dynamic_scalar_model(d3_data);
solve!(d3_model, cumulative_iteration_limit=0)
set_value!.(d3_model[:TAX][2010:2100], 0.1);
solve!(d3_model)
df3 = dyn_model_report(d3_model, d3_data; value_name = :Adjust_Capital)



###############
## Reporting ##
###############


df = vcat(
    df1,
    df2,
    df3
) 



p_i = df |>
    x -> subset(x, :variable => ByRow(==("invest"))) |>
    df -> plot(
        df,
        x=:time,
        y=:value,
        color=:model,
        Layout(
            title="Investment (% deviation from baseline)", 
            yaxis_title="% deviation from baseline", 
            xaxis_title="Year"
            ),
    )

p_o = df |>
    x -> subset(x, :variable => ByRow(==("output"))) |>
    df -> plot(
        df,
        x=:time,
        y=:value,
        color=:model,
        Layout(
            title="Output (% deviation from baseline)", 
            yaxis_title="% deviation from baseline", 
            xaxis_title="Year"
            ),
    )    


p_cons = df |>
    x -> subset(x, :variable => ByRow(==("cons"))) |>
    df -> plot(
        df,
        x=:time,
        y=:value,
        color=:model,
        Layout(
            title="Consumption (% deviation from baseline)", 
            yaxis_title="% deviation from baseline", 
            xaxis_title="Year"
            ),
    )    

    
p_capital = df |>
    x -> subset(x, :variable => ByRow(==("capital"))) |>
    df -> plot(
        df,
        x=:time,
        y=:value,
        color=:model,
        Layout(
            title="Capital (% deviation from baseline)", 
            yaxis_title="% deviation from baseline", 
            xaxis_title="Year"
            ),
    )






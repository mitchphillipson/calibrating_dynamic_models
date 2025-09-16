using MPSGE
using DataFrames
using PlotlyJS
import JuMP.Containers: DenseAxisArray



include("scalar_models.jl")

static_data = StaticData(
    200, # Y0
    180, # C0
    150, # LS0
    50,  # KS0
    20  # I0
    )

static_model = DYN_Static(static_data)
solve!(static_model, cumulative_iteration_limit=0)

set_value!(static_model[:TAX], 0.1);
solve!(static_model)


d1_data = DYN1_data(static_data);
d1_model = DYN_Scalar(d1_data);
solve!(d1_model, cumulative_iteration_limit=0)
set_value!.(d1_model[:TAX][2010:2100], 0.1);
solve!(d1_model)

df1 = dyn_model_report(d1_model, d1_data; value_name = :dyn1)



d2_data = DYN2_data(static_data);
d2_model = DYN_Scalar(d2_data);
solve!(d2_model, cumulative_iteration_limit=0)
set_value!.(d2_model[:TAX][2010:2100], 0.1);
solve!(d2_model)

df2 = dyn_model_report(d2_model, d2_data; value_name = :dyn2)



d3_data = DYN3_data(static_data);
d3_model = DYN_Scalar(d3_data);
solve!(d3_model, cumulative_iteration_limit=0)
set_value!.(d3_model[:TAX][2010:2100], 0.1);
solve!(d3_model)

df3 = dyn_model_report(d3_model, d3_data; value_name = :dyn3)



df = outerjoin(
    df1,
    df2,
    df3,
    on = [:time, :variable],
    makeunique=true,
) |>
x -> stack(x, Not(:time, :variable), variable_name = :model, value_name = :value)



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

p_i = df |>
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



df2 |>
    x -> subset(x, :variable => ByRow(==("output")))


value.(d2_model[:Y])




value(d2_model[:Y][2001])/d2_data.QREF[2001]-1

1.0203/d2_data.QREF[2001]-1



df1 |>
    x -> subset(x, :variable => ByRow(==("output")))

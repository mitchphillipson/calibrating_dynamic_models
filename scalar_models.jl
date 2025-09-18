"""
    dynamic_scalar_model(data::DynamicDataScalar)

Creates a dynamic scalar MPSGE model based on the provided `data`.
"""
function dynamic_scalar_model(data::DynamicDataScalar)

    Y0 = data.Y0
    C0 = data.C0
    LS0 = data.LS0
    KS0 = data.KS0
    I0 = data.I0

    time_periods = data.time_periods
    time_periods_horizon = data.time_periods_horizon

    G = data.G
    D = data.D
    K0 = data.K0
    rk0 = data.rk0

    QREF = data.QREF
    r = data.r
    pk0 = data.pk0
    PREF = data.PREF


    DYN1 = MPSGEModel()

    @parameter(DYN1, TAX[t=time_periods], 0)

    @sectors(DYN1, begin
        Y[t=time_periods], (start = QREF[t], description = "Production")
        I[t=time_periods], (start = QREF[t], description = "Investment")
        K[t=time_periods], (start = QREF[t], description = "Capital")
        C[t=time_periods], (start = QREF[t], description = "Consumption Index")
    end)

    @commodities(DYN1, begin
        PY[t=time_periods], (start = PREF[t], description = "Price index on output")
        PC[t=time_periods], (start = PREF[t], description = "Price index on consumption")
        RK[t=time_periods], (start = PREF[t], description = "Present value return to capital")
        PL[t=time_periods], (start = PREF[t], description = "Present value wage")
        PK[t=time_periods_horizon], (start = PREF[t]*pk0, description = "Price index on capital")
        #PKT -> PKT = PK[2101]
    end)

    @consumer(DYN1, RA)

    @auxiliary(DYN1, TCAP, start = I0*QREF[end] + K0*(1-D)*QREF[end])

    @production(DYN1, Y[t=time_periods], [t=0, s=1], begin
        @output(PY[t], Y0, t)
        @input(RK[t], KS0, s)
        @input(PL[t], LS0, s)
    end)

    @production(DYN1, I[t=time_periods], [t=0, s=0], begin
        @output(PK[t+1], I0, t)
        @input(PY[t], I0, s, taxes = [Tax(RA, TAX[t])])
    end)

    @production(DYN1, K[t=time_periods], [t=0, s=0], begin
        @output(PK[t+1], K0*(1-D), t)
        @output(RK[t], KS0, t)
        @input(PK[t], K0, s)
    end)

    @production(DYN1, C[t=time_periods], [t=0, s=0], begin
        @output(PC[t], C0, t)
        @input(PY[t], C0, s)
    end)

    @demand(DYN1, RA, begin
        @final_demand(PC[t=time_periods], C0*QREF[t], reference_price = PREF[t])
        @endowment(PK[time_periods[begin]], K0)
        @endowment(PL[t=time_periods], LS0*QREF[t])
        @endowment(PK[time_periods_horizon[end]], -TCAP)
    end)


    @aux_constraint(DYN1, TCAP, 
        C[time_periods[end]-1]*I[time_periods[end]] - I[time_periods[end]-1]*C[time_periods[end]]
    )

    return DYN1
end

"""
    dyn_model_report(model::MPSGEModel, data::DynamicDataScalar; value_name = :value)

Generates a DataFrame report for the dynamic model, showing percentage deviations 
from the baseline for key variables over time.

The DataFramw will contain four columns: `:time`, `:variable`, `:value`, and `:model`.

The `:variable` column indicates the economic variable (e.g., `:invest`, 
`:cons`, `:capital`, `:output`), the `:value` column contains the percentage
deviation from the baseline, and the `:model` column is a constant indicating
the name of the model (as specified by `value_name`).

Values of each variable are calculated as:

    100 * (value(model[:VAR][t]) / QREF[t] - 1)

where `VAR` is the variable of interest and `QREF[t]` is the baseline value as 
defined in `data.QREF`.
"""
function dyn_model_report(model::MPSGEModel, data::DynamicDataScalar; value_name = :value)
    time_periods = data.time_periods
    QREF = data.QREF

    df = DataFrame(
        time = time_periods,
        invest = [100*(value(model[:I][t])/QREF[t] - 1) for t in time_periods],
        cons = [100* ( value(model[:C][t])/QREF[t] - 1) for t in time_periods],
        capital = [100*(value(model[:K][t])/QREF[t] - 1) for t in time_periods],
        output = [100*(value(model[:Y][t])/QREF[t] - 1) for t in time_periods]
    ) |>
    x -> stack(x, Not(:time), variable_name = :variable, value_name = :value) |>
    x -> transform(x,
        :variable => ByRow(y-> value_name) => :model
    )

    return df
end


"""
    static_scalar_model(data::StaticData)

Creates a static scalar MPSGE model based on the provided `data`.
"""
function static_scalar_model(data::StaticData)

    Y0 = data.Y0
    C0 = data.C0
    LS0 = data.LS0
    KS0 = data.KS0
    I0 = data.I0

    SingleYear = MPSGEModel()

    @parameter(SingleYear, TAX, 0)

    @sectors(SingleYear, begin
        Y
        I
    end)

    @commodities(SingleYear, begin
        PY
        PL
        RK
        PK
    end)

    @consumer(SingleYear, RA)


    @production(SingleYear, Y, [t=0, s=1], begin
        @output(PY, Y0, t)
        @input(PL, LS0, s)
        @input(RK, KS0, s)
    end)

    @production(SingleYear, I, [t=0, s=0], begin
        @output(PK, I0, t)
        @input(PY, I0, s, taxes = [Tax(RA, TAX)])
    end)

    @demand(SingleYear, RA, begin
        @final_demand(PY, C0)
        @endowment(PL, LS0)
        @endowment(RK, KS0)
        @endowment(PK, -I0)
    end)


    return SingleYear
end
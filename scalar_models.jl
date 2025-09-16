using MPSGE
using DataFrames

import JuMP.Containers: DenseAxisArray

struct StaticData
    Y0::Float64
    C0::Float64
    LS0::Float64
    KS0::Float64
    I0::Float64
end


struct DYN_data
    Y0::Float64
    C0::Float64
    LS0::Float64
    KS0::Float64
    I0::Float64
    G::Float64
    D::Float64
    K0::Float64
    rk0::Float64
    r::Float64
    pk0::Float64

    time_periods::UnitRange{Int}
    time_periods_horizon::UnitRange{Int}
    QREF::DenseAxisArray{Float64, 1, Tuple{UnitRange{Int}}}
    PREF::DenseAxisArray{Float64, 1, Tuple{UnitRange{Int}}}
end


function DYN1_data(;
        start_year = 2000, 
        end_year = 2100, 
        Y0 = 200, 
        C0 = 180, 
        LS0 = 150, 
        KS0 = 50, 
        I0 = 20, 
        G = 0.02, 
        D = 0.02, 
        method = :DOC
    )
    time_periods = start_year:end_year
    time_periods_horizon = start_year:(end_year+1)

    K0 = I0/(G+D)

    rk0 = KS0/K0

    QREF = DenseAxisArray((1+G).^(eachindex(time_periods).-1), time_periods)
    r = rk0-D
    pk0 = 1+r
    PREF = DenseAxisArray((1 / (1+r)).^(eachindex(time_periods_horizon).-1), time_periods_horizon)
    #PREF = DenseAxisArray((1-r).^(eachindex(time_periods_horizon).-1), time_periods_horizon)

    ## As GAMS states ##
    if method == :GMS
        r = (rk0-D)/(1+rk0-D)
    end

    return DYN_data(Y0, C0, LS0, KS0, I0, G, D, K0, rk0, r, pk0, time_periods, time_periods_horizon, QREF, PREF)

end

function DYN1_data(
        static_data::StaticData; 
        start_year = 2000,
        end_year = 2100,
        G = 0.02,
        D = 0.02,
        method = :DOC
    )

    return DYN1_data(
        start_year = start_year,
        end_year = end_year,
        Y0 = static_data.Y0,
        C0 = static_data.C0,
        LS0 = static_data.LS0,
        KS0 = static_data.KS0,
        I0 = static_data.I0,
        G = G,
        D = D,
        method = method
    )

end


function DYN2_data(;
        start_year = 2000, 
        end_year = 2100, 
        Y0 = 200, 
        C0 = 180, 
        LS0 = 150, 
        KS0 = 50, 
        I0 = 20, 
        G = 0.02, 
        D = 0.02, 
        R = 0.05
    )
    time_periods = start_year:end_year
    time_periods_horizon = start_year:(end_year+1)

    pk0 = 1/(1-R)
    rk0 = (R-R*D+D)/(1-R)
    K0 = KS0/rk0
    adjusted_I0 = K0*(G+D)
    adjusted_C0 = C0 + I0 - adjusted_I0

    QREF = DenseAxisArray((1+G).^(eachindex(time_periods).-1), time_periods)
    PREF = DenseAxisArray((1-R).^(eachindex(time_periods_horizon).-1), time_periods_horizon)
    return DYN_data(Y0, adjusted_C0, LS0, KS0, adjusted_I0, G, D, K0, rk0, R, pk0, time_periods, time_periods_horizon, QREF, PREF)

end

function DYN2_data(
        static_data::StaticData; 
        start_year = 2000,
        end_year = 2100,
        G = 0.02,
        D = 0.02,
        R = 0.05
    )

    return DYN2_data(
        start_year = start_year,
        end_year = end_year,
        Y0 = static_data.Y0,
        C0 = static_data.C0,
        LS0 = static_data.LS0,
        KS0 = static_data.KS0,
        I0 = static_data.I0,
        G = G,
        D = D,
        R = R
    )

end


function DYN3_data(;
        start_year = 2000, 
        end_year = 2100, 
        Y0 = 200, 
        C0 = 180, 
        LS0 = 150, 
        KS0 = 50, 
        I0 = 20, 
        G = 0.02, 
        D = 0.02, 
        R = 0.05
    )
    time_periods = start_year:end_year
    time_periods_horizon = start_year:(end_year+1)

    pk0 = 1/(1-R)
    rk0 = (R-R*D+D)/(1-R)
    K0 = I0/(G+D)

    adjusted_KS0 = K0*rk0
    adjusted_LS0 = LS0 + KS0 - adjusted_KS0


    QREF = DenseAxisArray((1+G).^(eachindex(time_periods).-1), time_periods)
    PREF = DenseAxisArray((1-R).^(eachindex(time_periods_horizon).-1), time_periods_horizon)
    return DYN_data(Y0, C0, adjusted_LS0, adjusted_KS0, I0, G, D, K0, rk0, R, pk0, time_periods, time_periods_horizon, QREF, PREF)

end


function DYN3_data(
        static_data::StaticData; 
        start_year = 2000,
        end_year = 2100,
        G = 0.02,
        D = 0.02,
        R = 0.05
    )

    return DYN3_data(
        start_year = start_year,
        end_year = end_year,
        Y0 = static_data.Y0,
        C0 = static_data.C0,
        LS0 = static_data.LS0,
        KS0 = static_data.KS0,
        I0 = static_data.I0,
        G = G,
        D = D,
        R = R
    )

end


function DYN_Scalar(data::DYN_data)

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


function dyn_model_report(model::MPSGEModel, data::DYN_data; value_name = :value)
    time_periods = data.time_periods
    QREF = data.QREF

    df = DataFrame(
        time = time_periods,
        invest = [100*(value(model[:I][t])/QREF[t] - 1) for t in time_periods],
        cons = [100* ( value(model[:C][t])/QREF[t] - 1) for t in time_periods],
        capital = [100*(value(model[:K][t])/QREF[t] - 1) for t in time_periods],
        output = [100*(value(model[:Y][t])/QREF[t] - 1) for t in time_periods]
    ) |>
    x -> stack(x, Not(:time), variable_name = :variable, value_name = value_name)

    return df
end



function DYN_Static(data::StaticData)

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
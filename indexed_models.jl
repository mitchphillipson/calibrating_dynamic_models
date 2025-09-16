import JuMP
using Ipopt

import JuMP.Containers: DenseAxisArray

struct IndexedData
    Y0::Float64
    C0::Float64
    LS0::Float64
    KS0::Float64
    I0::Float64

    goods::Vector{Symbol}
    X::DenseAxisArray
    LD::DenseAxisArray
    KD::DenseAxisArray

    G::Float64
    R::Float64
    D::Float64
    K0::Float64
    rk0::Float64
    
    pk0::Float64

    time_periods::UnitRange{Int}
    time_periods_horizon::UnitRange{Int}
    QREF::DenseAxisArray{Float64, 1, Tuple{UnitRange{Int}}}
    PREF::DenseAxisArray{Float64, 1, Tuple{UnitRange{Int}}}

end


function IndexedData(;
        Y0 = 200,
        C0 = 180,
        LS0 = 150,
        KS0 = 50,
        I0 = 20,
        goods = [:A, :B],
        X = DenseAxisArray([120.0, 80.0], goods),
        LD = DenseAxisArray([100.0, 50.0], goods),
        KD = DenseAxisArray([20.0, 30.0], goods),
        G = 0.02,
        R = 0.05,
        D = 0.02,
        start_year = 2000,
        end_year = 2100,
    )


    time_periods = start_year:end_year
    time_periods_horizon = start_year:(end_year+1)
    K0 = I0/(G+D)
    rk0 = (R-R*D+D)/(1-R) 
    pk0 = 1/(1-R)

    
    QREF = DenseAxisArray((1+G).^(eachindex(time_periods).-1), time_periods)
    PREF = DenseAxisArray((1 / (1+R)).^(eachindex(time_periods_horizon).-1), time_periods_horizon)

    adjusted_KS0 = K0*rk0

    M = JuMP.Model(Ipopt.Optimizer)
    JuMP.@variables(M, begin
        VK[g=goods] >= 0 # Calibrated value of capital earnings
        VL[g=goods] >= 0 # Calibrated value of labor earnings
    end)

    JuMP.@objective(M, Min, sum(1/KD[g]*(VK[g] - KD[g])^2 for g∈goods))

    JuMP.@constraints(M, begin
        VABAL[g=goods], VK[g] + VL[g] == KD[g] + LD[g]
        VKBAL, sum(VK[g] for g∈goods) == adjusted_KS0
    end)

    JuMP.optimize!(M)

    adjusted_LD = JuMP.value.(VL)
    adjusted_KD = JuMP.value.(VK)


    return IndexedData(Y0, C0, LS0, adjusted_KS0, I0, goods, X, adjusted_LD, adjusted_KD, G, R, D, K0, rk0, pk0, time_periods, time_periods_horizon, QREF, PREF)
end

function IndexedData(
    static_data::StaticData;
    goods = [:A, :B],
    X = DenseAxisArray([120.0, 80.0], goods),
    LD = DenseAxisArray([100.0, 50.0], goods),
    KD = DenseAxisArray([20.0, 30.0], goods),
    G = 0.02,
    R = 0.05,
    D = 0.02,
    start_year = 2000,
    end_year = 2100,
    )

    return IndexedData(
        static_data.Y0,
        static_data.C0,
        static_data.LS0,
        static_data.KS0,
        static_data.I0,
        goods,
        X,
        LD,
        KD,
        G,
        R,
        D,
        start_year = start_year,
        end_year = end_year
    )
end



function DYN_Indexed(data::IndexedData)

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
    r = data.R
    pk0 = data.pk0
    PREF = data.PREF

    goods = data.goods

    X0 = data.X
    LD = data.LD
    KD = data.KD


    DYN = MPSGEModel()

    @parameter(DYN, TAX[t=time_periods], 0)

    @sectors(DYN, begin
        Y[t=time_periods], (start = QREF[t], description = "Macro Output (transitory utility)")
        X[g=goods, t=time_periods], (start = QREF[t], description = "Production")
        I[t=time_periods], (start = QREF[t], description = "Investment")
        K[t=time_periods], (start = QREF[t], description = "Capital")
        C[t=time_periods], (start = QREF[t], description = "Consumption Index")
    end)

    @commodities(DYN, begin
        PY[t=time_periods], (start = PREF[t], description = "Price index on output")
        PX[g=goods, t=time_periods], (start = PREF[t], description = "Price index on sector output")
        PC[t=time_periods], (start = PREF[t], description = "Price index on consumption")
        RK[t=time_periods], (start = PREF[t], description = "Present value return to capital")
        PL[t=time_periods], (start = PREF[t], description = "Present value wage")
        PK[t=time_periods_horizon], (start = PREF[t]*pk0, description = "Price index on capital")
        #PKT -> PKT = PK[2101]
    end)

    @consumer(DYN, RA)

    @auxiliary(DYN, TCAP, start = I0*QREF[end] + K0*(1-D)*QREF[end])

    @production(DYN, Y[t=time_periods], [t=0, s=1], begin
        @output(PY[t], Y0, t)
        @input(PX[g=goods, t], X0[g], s)
    end)


    @production(DYN, X[g=goods, t=time_periods], [t=0, s=1], begin
        @output(PX[g, t], X0[g], t)
        @input(RK[t], KD[g], s)
        @input(PL[t], LD[g], s)
    end)

    @production(DYN, I[t=time_periods], [t=0, s=0], begin
        @output(PK[t+1], I0, t)
        @input(PY[t], I0, s, taxes = [Tax(RA, TAX[t])])
    end)

    @production(DYN, K[t=time_periods], [t=0, s=0], begin
        @output(PK[t+1], K0*(1-D), t)
        @output(RK[t], KS0, t)
        @input(PK[t], K0, s)
    end)

    @production(DYN, C[t=time_periods], [t=0, s=0], begin
        @output(PC[t], C0, t)
        @input(PY[t], C0, s)
    end)

    @demand(DYN, RA, begin
        @final_demand(PC[t=time_periods], C0*QREF[t], reference_price = PREF[t])
        @endowment(PK[time_periods[begin]], K0)
        @endowment(PL[t=time_periods], sum(LD[g] for g in goods)*QREF[t])
        @endowment(PK[time_periods_horizon[end]], -TCAP)
    end)


    @aux_constraint(DYN, TCAP, 
        C[time_periods[end-1]]*I[time_periods[end]] - I[time_periods[end-1]]*C[time_periods[end]]
    )

    return DYN
end



d4_data = IndexedData();
d4_model = DYN_Indexed(d4_data);
solve!(d4_model, cumulative_iteration_limit=0)

generate_report(d4_model) |>
    x -> sort(x, :margin)


d4_data.K0
d4_data.KS0
d4_data.D

start_value.(d4_model[:RK])


production(d4_model[:K][2000])


value(d4_model[:PK][2001])
value(d4_model[:PK][2000])

value(d4_model[:RK][2000])


zero_profit(d4_model[:K][2000])


value(d4_model[:PK][2000])*d4_data.K0

d4_data.KS0*value(d4_model[:RK][2000])
d4_data.K0*(1-d4_data.D)*value(d4_model[:PK][2000])


PK = d4_model[:PK]
RK = d4_model[:RK]
K0 = d4_data.K0
KS0 = d4_data.KS0
D = d4_data.D

value(500*PK[2000])

value(.931*PK[2001] + .096*RK[2000])

526.31577*.931
tot = K0*(1-D) + KS0

K0*(1-D)*value(PK[2001])+ KS0*value(RK[2000])

K0*(1-D) + KS0

start_value(PK[2000])^2*K0

D
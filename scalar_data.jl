"""
    StaticData

A structure to hold static model data. Fields:

    - `Y0`  : Initial output
    - `C0`  : Initial consumption
    - `LS0` : Initial labor supply
    - `KS0` : Initial capital stock
    - `I0`  : Initial investment
"""
struct StaticData
    Y0::Float64
    C0::Float64
    LS0::Float64
    KS0::Float64
    I0::Float64
end


"""
    DynamicDataScalar

A structure to hold dynamic model data. Fields:

    - `Y0` : Initial output
    - `C0` : Initial consumption
    - `LS0` : Initial labor supply
    - `KS0` : Initial capital stock
    - `I0` : Initial investment
    - `G` : Growth rate of technology
    - `D` : Depreciation rate of capital
    - `K0` : Initial capital (efficiency units)
    - `rk0` : Initial rental rate of capital (efficiency units)
    - `r` : Interest rate
    - `pk0` : Initial price of capital
    - `time_periods` : Time periods for the model (e.g., years 2000 to 2100)
    - `time_periods_horizon` : Time periods including the horizon (e.g., years 2000 to 2101)
    - `QREF` : Reference output levels for each time period
    - `PREF` : Reference price levels for each time period in the horizon

It's recommended to use the provided functions `DYN1_data`, `dynamic_data_2`, or 
`dynamic_data_3` to create instances of this structure with appropriate parameters.
"""
struct DynamicDataScalar
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

"""
    dynamic_data_1(;
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

    dynamic_data_1(
        static_data::StaticData; 
        start_year = 2000,
        end_year = 2100,
        G = 0.02,
        D = 0.02,
        method = :DOC
    )

Create data need for the first dynamic model calibration approach. In this approach,
the interest rate is derived from the rental rate of capital and depreciation rate.

``
K0 = \\frac{I0}{G+D}
``

``
r = \\frac{KS0}{K0} - D
``
"""
function dynamic_data_1(;
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

    return DynamicDataScalar(Y0, C0, LS0, KS0, I0, G, D, K0, rk0, r, pk0, time_periods, time_periods_horizon, QREF, PREF)

end

function dynamic_data_1(
        static_data::StaticData; 
        start_year = 2000,
        end_year = 2100,
        G = 0.02,
        D = 0.02,
        method = :DOC
    )

    return dynamic_data_1(
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


function dynamic_data_2_gams(;
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
    return DynamicDataScalar(Y0, adjusted_C0, LS0, KS0, adjusted_I0, G, D, K0, rk0, R, pk0, time_periods, time_periods_horizon, QREF, PREF)

end



function dynamic_data_2_gams(
        static_data::StaticData; 
        start_year = 2000,
        end_year = 2100,
        G = 0.02,
        D = 0.02,
        R = 0.05
    )

    return dynamic_data_2_gams(
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

"""
    dynamic_data_2(;
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

    dynamic_data_2(
        static_data::StaticData; 
        start_year = 2000,
        end_year = 2100,
        G = 0.02,
        D = 0.02,
        R = 0.05
    )

Create data need for the second dynamic model calibration approach. In this approach,
the interest rate is given and we adjust investment

``
K0' = \\frac{KS0}{R+D}
``

``
I0' = K0'(G+D)
``

``
C0' = C0 + I0 - I0'
``
"""
function dynamic_data_2(;
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

    pk0 = 1+R
    rk0 = R+D #(R-R*D+D)/(1-R)
    K0 = KS0/rk0
    adjusted_I0 = K0*(G+D)
    adjusted_C0 = C0 + I0 - adjusted_I0

    QREF = DenseAxisArray((1+G).^(eachindex(time_periods).-1), time_periods)
    PREF = DenseAxisArray((1/(1+R)).^(eachindex(time_periods_horizon).-1), time_periods_horizon)
    return DynamicDataScalar(Y0, adjusted_C0, LS0, KS0, adjusted_I0, G, D, K0, rk0, R, pk0, time_periods, time_periods_horizon, QREF, PREF)

end

function dynamic_data_2(
        static_data::StaticData; 
        start_year = 2000,
        end_year = 2100,
        G = 0.02,
        D = 0.02,
        R = 0.05
    )

    return dynamic_data_2(
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



"""
    dynamic_data_3(;
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

    dynamic_data_3(
        static_data::StaticData; 
        start_year = 2000,
        end_year = 2100,
        G = 0.02,
        D = 0.02,
        R = 0.05
    )

Create data need for the third dynamic model calibration approach. In this approach,
the interest rate is given and we adjust investment

``
KS0' = K0(R+D)
``

``
LS0' = LS0 + KS0 - KS0'
``
"""
function dynamic_data_3(;
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

    pk0 = 1+R
    rk0 = R+D #(R-R*D+D)/(1-R)
    K0 = I0/(G+D)

    adjusted_KS0 = K0*rk0
    adjusted_LS0 = LS0 + KS0 - adjusted_KS0


    QREF = DenseAxisArray((1+G).^(eachindex(time_periods).-1), time_periods)
    PREF = DenseAxisArray((1-R).^(eachindex(time_periods_horizon).-1), time_periods_horizon)
    return DynamicDataScalar(Y0, C0, adjusted_LS0, adjusted_KS0, I0, G, D, K0, rk0, R, pk0, time_periods, time_periods_horizon, QREF, PREF)

end


function dynamic_data_3(
        static_data::StaticData; 
        start_year = 2000,
        end_year = 2100,
        G = 0.02,
        D = 0.02,
        R = 0.05
    )

    return dynamic_data_3(
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

function dynamic_data_3_gams(;
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
    return DynamicDataScalar(Y0, C0, adjusted_LS0, adjusted_KS0, I0, G, D, K0, rk0, R, pk0, time_periods, time_periods_horizon, QREF, PREF)

end


function dynamic_data_3_gams(
        static_data::StaticData; 
        start_year = 2000,
        end_year = 2100,
        G = 0.02,
        D = 0.02,
        R = 0.05
    )

    return dynamic_data_3_gams(
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


#### Cache structures
abstract type AbstractWeissingerResults end
struct WeissingerResults{T} <: AbstractWeissingerResults
    # Preallocate
    Nelempan::Vector{Int64}
    ielem::Vector{Int64}
    AIC::Matrix{T}
    ipiv::Vector{Int64}
    gamma::Matrix{T}
    w::Matrix{T}
    cl::Matrix{T}
    clcbar::Vector{T}
    CL::Vector{T}
    Cl::Vector{T}
    Cm::Vector{T}
    CDff::Vector{T}
    CDp::Vector{Vector{T}}
    CDi::Vector{Vector{T}}
    CD::Vector{T}

    # Grid arrays
    xc ::Vector{T}
    yc ::Vector{T}
    xr ::Vector{T}
    yr ::Vector{T}
    xt ::Vector{T}
    yt ::Vector{T}
    xctl ::Vector{T}
    xrc ::Vector{T}
    xtc ::Vector{T}
    chord ::Vector{T}
    paneli::Vector{T}
    area::Vector{T}
    function WeissingerResults(T, Nelempan, Ncontrols=length(Nelempan))
        ielem        = cumsum(Nelempan) - Nelempan # starting index of each element
        outlen       = length(deriv) + Ncontrols
        Nelem        = length(Nelempan)
        Npanhalfwing = sum(Nelempan)
        Npantotal    = 2 * Npanhalfwing # include mirrored side

        # Vortex solution arrays
        AIC    = zeros(T, Npantotal,Npantotal)
        ipiv   = zeros(Int64, Npantotal)
        gamma  = zeros(T, Npantotal, outlen)
        w      = zeros(T, Npantotal, outlen)
        cl     = zeros(T, Npantotal, outlen)
        clcbar = zeros(T, Npantotal)

        # Outputs
        CL      = zeros(T, outlen)
        Cl      = zeros(T, outlen)
        Cm      = zeros(T, outlen)
        CDff    = zeros(T, outlen)
        CD      = zeros(T, outlen)
        CDp     = [zeros(T, outlen) for _=1:Nelem]
        CDi     = [zeros(T, outlen) for _=1:Nelem]
    
        # Grid arrays
        xc     = zeros(T, Npanhalfwing)
        yc     = zeros(T, Npanhalfwing)
        xr     = zeros(T, Npanhalfwing)
        yr     = zeros(T, Npanhalfwing)
        xt     = zeros(T, Npanhalfwing)
        yt     = zeros(T, Npanhalfwing)
        xctl   = zeros(T, Npanhalfwing)
        xrc    = zeros(T, Npanhalfwing)
        xtc    = zeros(T, Npanhalfwing)
        chord  = zeros(T, Npanhalfwing)
        paneli = zeros(T, Npanhalfwing)
        area   = zeros(T, Nelem)
        new{T}( Nelempan, ielem, AIC, ipiv, gamma, w, cl, clcbar, CL,
                    Cl, Cm, CDff, CDp, CDi, CD, xc, yc, xr, yr, xt, yt, 
                    xctl, xrc, xtc, chord, paneli, area)
    end
end

struct WeissingerResultsAD{T,TF} <: AbstractWeissingerResults
    primal::WeissingerResults{T}
    dual::WeissingerResults{TF}
    Nelempan::Vector{Int64}
    ipiv::Vector{Int64}
    function WeissingerResultsAD(T, TF, Nelempan)
        new{T, TF}(WeissingerResults(T, Nelempan), WeissingerResults(TF, Nelempan),Nelempan,zeros(Int64,2*sum(Nelempan)))
    end
end
@inline get_array(wres::WeissingerResultsAD, key::Symbol, ::Type{ForwardDiff.Dual{Tg,V,N}} where {Tg,V,N}) = getfield(wres.dual, key)
@inline get_array(wres::WeissingerResultsAD, key::Symbol, ::Type{<:Number}) = getfield(wres.primal, key)
@inline get_array(wres::WeissingerResults, key::Symbol, ::Type{<:Number}) = getfield(wres, key)
@inline get_results(wres::WeissingerResultsAD, ::Type{ForwardDiff.Dual{Tg,V,N}} where {Tg,V,N}) = wres.dual
@inline get_results(wres::WeissingerResultsAD, ::Type{<:Number}) = wres.primal
@inline get_results(wres::WeissingerResults, ::Type{<:Number}) = wres


function primalize(result::WeissingerResults, resultAD::WeissingerResultsAD)
    for k = propertynames(resultAD.dual)
        v = getfield(resultAD.dual,k)
        if eltype(v) <: ForwardDiff.Dual
            vp = getfield(result,k)
            vp .= ForwardDiff.value.(v)
        end
    end
    return result
end
primalize(resultAD::WeissingerResultsAD) = primalize(resultAD.primal, resultAD)

function Base.copy(wres::WeissingerResults{T}) where T
    wnew = WeissingerResults(T, wres.Nelempan)
    
    wnew.AIC    .= wres.AIC
    wnew.ipiv   .= wres.ipiv
    wnew.gamma  .= wres.gamma
    wnew.w      .= wres.w
    wnew.cl     .= wres.cl
    wnew.clcbar .= wres.clcbar
    wnew.CL     .= wres.CL
    wnew.Cl     .= wres.Cl
    wnew.Cm     .= wres.Cm
    wnew.CDff   .= wres.CDff
    wnew.CDp    .= wres.CDp
    wnew.CDi    .= wres.CDi
    wnew.CD     .= wres.CD

    wnew.xc     .= wres.xc
    wnew.yc     .= wres.yc
    wnew.xr     .= wres.xr
    wnew.yr     .= wres.yr
    wnew.xt     .= wres.xt
    wnew.yt     .= wres.yt
    wnew.xctl   .= wres.xctl
    wnew.xrc    .= wres.xrc
    wnew.xtc    .= wres.xtc
    wnew.chord  .= wres.chord
    wnew.paneli .= wres.paneli
    wnew.area   .= wres.area
    return wnew
end
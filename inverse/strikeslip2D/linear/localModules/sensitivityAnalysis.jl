# Raw ∂J/∂p, logarithmic p·∂J/∂p and relative (p/J)·∂J/∂p sensitivities of η, ρ, G, K
#
# - Gradients are derivatives w.r.t. the discrete field values.
# - The adjoint of a step holds the earlier states fixed.
# - log10η = ∂J/∂log10(η) = ln(10)·η·∂J/∂η: linearised change in J per decade of viscosity.
using Statistics

const SENSITIVITY_LABELS = (η=L"\eta",log10η=L"\log_{10}\eta",ρ=L"\rho",G=L"G",K=L"K")

function checked_field(A, coords, name)
    B = Array(A)
    ndims(B)==2 || error("$name must be 2D; got size $(size(B))")
    expected = (length(coords[1]),length(coords[2]))
    size(B)==expected || error("$name: size $(size(B)) != grid size $expected")
    return B
end

function summed_phase_field(A, coords, name)
    B = Array(A)
    ndims(B)==3 || error("$name must have (phase,x,y) dimensions; got $(size(B))")
    return checked_field(dropdims(sum(B;dims=1);dims=1),coords,name)
end

# J from the forward solution: Σ w·Vx for a weighted observation, else Σ Vx over the box
function evaluate_observation(stokes,grid,observation)
    observation.field==:Vx || error("Only observation.field=:Vx is implemented")
    if haskey(observation,:weights)
        terms = Array(observation.weights) .* Array(stokes.V.Vx)
        J = sum(terms); Jabs = sum(abs,terms)
        cancellation = iszero(Jabs) ? 0.0 : abs(J)/Jabs
        return (J=J,Jabs=Jabs,cancellation=cancellation,npoints=count(!iszero,observation.weights))
    end
    coords = (grid.xvi[1],grid.xci[2])
    Vx = checked_field(stokes.V.Vx[:,2:end-1],coords,"Vx")
    xc,yc = observation.center
    hx,hy = observation.half_width
    ix = findall(x->abs(x-xc)<=hx,coords[1])
    iy = findall(y->abs(y-yc)<=hy,coords[2])
    isempty(ix) && error("Observation box has no Vx x-nodes")
    isempty(iy) && error("Observation box has no Vx y-nodes")
    observed = Vx[ix,iy]
    J = sum(observed)
    Jabs = sum(abs,observed)
    cancellation = iszero(Jabs) ? 0.0 : abs(J)/Jabs
    return (J=J,Jabs=Jabs,cancellation=cancellation,npoints=length(observed))
end

# η on the vertices; ρ, G, K on the centres (center gradients include the folded vertex part)
function collect_sensitivity_fields(stokes,stokes_ad,gradients,grid;
                                    ρ_matrix,ρ_weak,G,K)
    ρ_matrix==ρ_weak || error("Density field construction currently assumes equal densities")
    xci,xvi = grid.xci,grid.xvi
    ηv = checked_field(stokes.viscosity.ηv,xvi,"forward ηv")
    raw = (
        η = checked_field(stokes_ad.viscosity.ηv,xvi,"dJ/dηv"),
        ρ = checked_field(stokes_ad.ρ,xci,"dJ/dρ"),
        G = summed_phase_field(gradients.G.center,xci,"dJ/dG center"),
        K = summed_phase_field(gradients.Kb.center,xci,"dJ/dKb center"),
    )
    params = (η=ηv,ρ=fill(Float64(ρ_matrix),size(raw.ρ)),
              G=fill(Float64(G),size(raw.G)),K=fill(isfinite(K) ? Float64(K) : 0.0,size(raw.K)))  # K = Inf (incompressible): K ∂J/∂K := 0

    # p ∂J/∂p = ∂J/∂ln(p): change in J per relative (e-fold) change in p
    log_e = NamedTuple{keys(raw)}(map(k->raw[k].*params[k],keys(raw)))
    logarithmic = (η=log_e.η, log10η=log(10).*log_e.η,
                   ρ=log_e.ρ, G=log_e.G, K=log_e.K)
    coords = (η=xvi,log10η=xvi,ρ=xci,G=xci,K=xci)
    return (raw=raw,logarithmic=logarithmic,params=params,coords=coords)
end

# (p/J) ∂J/∂p; skipped when J cancels (|J| ≪ Σ|terms|)
function relative_fields(logarithmic,objective;min_cancellation=1e-3)
    if !isfinite(objective.J) || objective.Jabs==0 || objective.cancellation<min_cancellation
        @warn "Relative sensitivities skipped: objective cancels or is zero" J=objective.J Jabs=objective.Jabs cancellation=objective.cancellation
        return nothing
    end
    return NamedTuple{keys(logarithmic)}(map(k->logarithmic[k]./objective.J,keys(logarithmic)))
end

function field_stats(fields)
    return NamedTuple{keys(fields)}(map(k->begin
        a=filter(isfinite,vec(fields[k]))
        isempty(a) ? (min=NaN,max=NaN,maxabs=NaN,rms=NaN,sum=NaN) :
        (min=minimum(a),max=maximum(a),maxabs=maximum(abs,a),
         rms=sqrt(mean(abs2,a)),sum=sum(a))
    end,keys(fields)))
end

# =============================================================================
# SENSITIVITY MAPS
# =============================================================================

const MAP_PX = 330   # panel width in px; height follows the domain aspect

# x rounded up to 1, 1.2, 1.5, 2, 2.5, 3, 4, 5, 6 or 8 × 10^e
function nice_limit(x)
    (isfinite(x) && x > 0) || return 1.0
    e = floor(log10(x))
    m = x/10.0^e
    for s in (1.0,1.2,1.5,2.0,2.5,3.0,4.0,5.0,6.0,8.0,10.0)
        m <= s*(1+1e-9) && return s*10.0^e
    end
    return 10.0^(e+1)
end

# Five ticks in ±vmax without the common power of ten, which goes into the colorbar label
function symmetric_ticks(vmax)
    e = floor(Int,log10(vmax))
    vals = [-vmax,-vmax/2,0.0,vmax/2,vmax]
    labels = [v==0 ? "0" : @sprintf("%.3g",v/10.0^e) for v in vals]
    return (vals,labels),e
end

function colorbar_label(e;quantity="",unit="")
    s = isempty(quantity) ? "" : quantity*"\\quad "
    s *= "\\times 10^{$e}"
    isempty(unit) || (s *= "\\;"*unit)
    return latexstring(s)
end

# Limits over vertex and centre grids together, so the panels align
function common_limits(coords,ks)
    xl = (minimum(first(coords[p][1]) for p in ks),maximum(last(coords[p][1]) for p in ks))
    yl = (minimum(first(coords[p][2]) for p in ks),maximum(last(coords[p][2]) for p in ks))
    return xl,yl
end

function sensitivity_map_axis(pos,i,p,xl,yl)
    hpx = round(Int,MAP_PX*(yl[2]-yl[1])/(xl[2]-xl[1]))
    ax = Axis(pos;title=SENSITIVITY_LABELS[p],titlesize=24,
              xlabel="x [m]",ylabel="y [m]",
              width=MAP_PX,height=hpx,
              xgridvisible=false,ygridvisible=false,
              xticklabelsize=15,yticklabelsize=15)
    i>1 && hideydecorations!(ax;ticks=false,minorticks=false)
    return ax
end

function decorate_map!(ax,xl,yl;x_min,x_max,y_min,y_max,
                       observation_center,observation_half_width)
    add_shear_zone!(ax,x_min,x_max,y_min,y_max)
    add_observation_box!(ax,observation_center,observation_half_width)
    limits!(ax,xl...,yl...)
    return nothing
end

# One colour scale per panel
function plot_scaled_sensitivity_maps(fields,coords,figdir;
                                      name,x_min,x_max,y_min,y_max,
                                      observation_center,observation_half_width,
                                      unit="")
    fields===nothing && return nothing
    ks = collect(keys(fields))
    xl,yl = common_limits(coords,ks)
    deco = (x_min=x_min,x_max=x_max,y_min=y_min,y_max=y_max,
            observation_center=observation_center,
            observation_half_width=observation_half_width)

    fig = Figure(fontsize=18,figure_padding=(10,25,10,10))
    for (i,p) in enumerate(ks)
        A = fields[p]
        x,y = coords[p]
        amax = maximum(abs,filter(isfinite,vec(A));init=0.0)
        zero_field = iszero(amax)
        vmax = zero_field ? 1.0 : nice_limit(amax)

        gl = fig[1,i] = GridLayout()
        ax = sensitivity_map_axis(gl[1,1],i,p,xl,yl)
        hm = heatmap!(ax,x,y,A;colormap=:balance,colorrange=(-vmax,vmax))
        decorate_map!(ax,xl,yl;deco...)
        zero_field && text!(ax,0.5,0.5;text=L"\equiv 0",space=:relative,
                            align=(:center,:center),fontsize=30,color=:gray40)

        ticks,e = symmetric_ticks(vmax)
        Colorbar(gl[2,1],hm;vertical=false,flipaxis=false,
                 width=Relative(0.9),height=10,
                 ticks=ticks,ticklabelsize=14,labelsize=17,
                 label=zero_field ? "" : colorbar_label(e;unit=unit))
        rowgap!(gl,8)
    end

    colgap!(fig.layout,22)
    resize_to_layout!(fig)

    save_figure(figdir,name,fig)
    return fig
end

function plot_maximum_sensitivities(relative,figdir;step=nothing)
    relative===nothing && return nothing
    stats=field_stats(relative)
    names=collect(keys(relative))
    n=length(names)
    vals=[stats[k].maxabs for k in names]
    fig=Figure(size=(200n+150,600),fontsize=22)
    ax=Axis(fig[1,1];ylabel=L"\max\,|(p/J)\,\partial J/\partial p|",yscale=log10,yticks=decade_ticks,ytickformat=decade_labels,
            xticks=(1:n,[SENSITIVITY_LABELS[k] for k in names]))
    positive=findall(v->isfinite(v)&&v>0,vals)
    if !isempty(positive)
        barplot!(ax,positive,vals[positive])
    else
        @warn "No positive relative sensitivities to plot"
    end
    save_figure(figdir,"sensitivity_max_relative"*step_tag(step),fig)
    return fig
end

function run_sensitivity_analysis(stokes,stokes_ad,gradients,grid,observation,figdir;
                                  ρ_matrix,ρ_weak,G,K,x_min,x_max,y_min,y_max,
                                  min_cancellation=1e-3,step=nothing,t=nothing,τM=nothing,plot=(:logarithmic,:relative),combined=false)
    mkpath(figdir)
    tag = step_tag(step)
    objective=evaluate_observation(stokes,grid,observation)
    log_at(2) && !isnothing(step) && @printf("\n---- Sensitivity analysis, step %d ----",step)
    log_at(2) && @printf("\nSensitivity objective J=% .8e, sum|Vx|=% .8e, cancellation=%.3e\n",
            objective.J,objective.Jabs,objective.cancellation)
    data=collect_sensitivity_fields(stokes,stokes_ad,gradients,grid;
                                    ρ_matrix,ρ_weak,G,K)
    relative=relative_fields(data.logarithmic,objective;min_cancellation)
    common=(x_min=x_min,x_max=x_max,y_min=y_min,y_max=y_max,
            observation_center=haskey(observation,:center) ? observation.center : nothing,
            observation_half_width=haskey(observation,:half_width) ? observation.half_width : nothing)

    combined && plot_sensitivities_combined(stokes_ad,gradients,grid,data.logarithmic,data.coords,figdir;
                                            step=step,t=t,τM=τM,common...)

    if :logarithmic in plot
        plot_scaled_sensitivity_maps(data.logarithmic,data.coords,figdir;
                                    name="sensitivity_logarithmic"*tag,unit=raw"\mathrm{m/s}",common...)
        plot_logarithmic_sensitivities_shared(data.logarithmic,data.coords,figdir;step=step,common...)
    end
    if :relative in plot
        plot_scaled_sensitivity_maps(relative,data.coords,figdir;
                                    name="sensitivity_relative"*tag,common...)
        plot_relative_sensitivities_shared(relative,data.coords,figdir;step=step,common...)
        plot_maximum_sensitivities(relative,figdir;step=step)
    end

    return (objective=objective,raw=data.raw,logarithmic=data.logarithmic,
            relative=relative,coords=data.coords)
end

# One colour scale for all panels
function plot_shared_sensitivity_maps(fields,coords,figdir;
                                      name,quantity,unit="",
                                      x_min,x_max,y_min,y_max,
                                      observation_center,observation_half_width)
    fields===nothing && return nothing
    ks = collect(keys(fields))
    n = length(ks)
    xl,yl = common_limits(coords,ks)
    deco = (x_min=x_min,x_max=x_max,y_min=y_min,y_max=y_max,
            observation_center=observation_center,
            observation_half_width=observation_half_width)

    Smax = maximum(maximum(abs,filter(isfinite,vec(fields[p]));init=0.0) for p in ks)
    vmax = iszero(Smax) ? 1.0 : nice_limit(Smax)
    log_at(2) && @printf("\n%s shared color scale: ±%.6e (data max %.6e)\n",name,vmax,Smax)

    fig = Figure(fontsize=18,figure_padding=(10,25,10,10))
    hm = nothing
    for (i,p) in enumerate(ks)
        A = fields[p]
        x,y = coords[p]
        ax = sensitivity_map_axis(fig[1,i],i,p,xl,yl)
        hm = heatmap!(ax,x,y,A;colormap=:balance,colorrange=(-vmax,vmax))
        decorate_map!(ax,xl,yl;deco...)
    end

    ticks,e = symmetric_ticks(vmax)
    Colorbar(fig[2,1:n],hm;vertical=false,flipaxis=false,
             width=Relative(0.45),height=12,
             ticks=ticks,ticklabelsize=15,labelsize=19,
             label=colorbar_label(e;quantity=quantity,unit=unit))

    colgap!(fig.layout,22)
    rowgap!(fig.layout,1,12)
    resize_to_layout!(fig)

    save_figure(figdir,name,fig)
    return fig
end

function plot_relative_sensitivities_shared(fields,coords,figdir;step=nothing,kwargs...)
    plot_shared_sensitivity_maps(fields,coords,figdir;
        name="sensitivity_relative_shared"*step_tag(step),
        quantity=raw"(p/J)\,\partial J/\partial p",
        kwargs...)
end

function plot_logarithmic_sensitivities_shared(fields,coords,figdir;step=nothing,kwargs...)
    plot_shared_sensitivity_maps(fields,coords,figdir;
        name="sensitivity_logarithmic_shared"*step_tag(step),
        quantity=raw"p\,\partial J/\partial p",
        unit=raw"\mathrm{m/s}",
        kwargs...)
end
# Figures of the forward fields and the sensitivities: functions copied from
# ../simpleShear/common/plotting.jl; only sensitivity_fields plots ∂J/∂G per cell (center)

# FIG_SETTINGS.show[]: also display every saved figure (slow in long loops)
const FIG_SETTINGS = (show=Ref(true),)
step_tag(step) = isnothing(step) ? "" : @sprintf("_step%04d",step)

# Frames (names ending in _stepNNNN) go into one subfolder per figure type;
# figures without a step tag stay in figdir.
const FRAME_FOLDERS = Dict(
    # forward
    "forward"=>"fields","velocity_field"=>"velocity",
    "slip_rate"=>"slip_rate","transient"=>"transient",
    # adjoint
    "adjoint"=>"adjoint_fields","sensitivity_heatmap"=>"dJdp_heatmap",
    "sensitivity_contours"=>"dJdp_contour","sensitivity_combined"=>"combined","sensitivity_creep"=>"creep","sensitivity_all"=>"all","sensitivity_all_heatmap"=>"all_heatmap",
    "sensitivity_2x2"=>"2x2","sensitivity_2x2_heatmap"=>"2x2_heatmap","sensitivity_logarithmic"=>"logarithmic",
    "sensitivity_logarithmic_shared"=>"logarithmic_shared","sensitivity_relative"=>"relative",
    "sensitivity_relative_shared"=>"relative_shared","sensitivity_max_relative"=>"max_relative",
    # viscosity inversion
    "viscosity_inversion"=>"viscosity_inversion")
const STEP_SUFFIX = r"_step\d+$"

function frame_dir(figdir,name)
    occursin(STEP_SUFFIX,name) || return figdir
    base = replace(name,STEP_SUFFIX=>"")
    return joinpath(figdir,get(FRAME_FOLDERS,base,base))
end

# The model folder lives on iCloud Drive, where closing a file occasionally times out
# (SystemError: close: Operation timed out); retry instead of losing the whole run.
function save_retry(path,fig;attempts=5,kwargs...)
    for k in 1:attempts
        try
            return save(path,fig;kwargs...)
        catch err
            (err isa SystemError && k < attempts) || rethrow()
            @warn "Saving $(basename(path)) failed ($(err.prefix)), retrying" attempt=k
            sleep(2.0k)
        end
    end
end

function save_figure(figdir,name,fig;px_per_unit=2)
    figdir = frame_dir(figdir,name)
    mkpath(figdir)
    save_retry(joinpath(figdir,"$name.png"),fig;px_per_unit=px_per_unit)
    FIG_SETTINGS.show[] && display(fig)
    return fig
end

function crange(A; diverging=false)
    a = filter(isfinite,vec(Array(A)))
    isempty(a) && return (-1.0,1.0)
    if diverging
        m = maximum(abs,a)
        return iszero(m) ? (-1.0,1.0) : (-m,m)
    end
    lo,hi = extrema(a)
    δ = iszero(lo) ? 1.0 : 0.01abs(lo)
    return lo == hi ? (lo-δ,hi+δ) : (lo,hi)
end

function panel!(fig,row,col,title,x,A; colormap=:batlow,diverging=false,nrows=2)
    ax = Axis(fig[row,2col-1];aspect=DataAspect(),title=title,
              xlabel=row==nrows ? "x" : "",ylabel=col==1 ? "y" : "",
              titlesize=20,xlabelsize=18,ylabelsize=18,xticklabelsize=15,yticklabelsize=15)
    hm = heatmap!(ax,x...,Array(A);colormap=colormap,colorrange=crange(A;diverging=diverging))
    Colorbar(fig[row,2col],hm;width=18,ticklabelsize=14)
    row<nrows && hidexdecorations!(ax;grid=false)
    col>1 && hideydecorations!(ax;grid=false)
    return ax
end

function add_shear_zone!(ax,x_min,x_max,y_min,y_max)
    lines!(ax,[x_min,x_max,x_max,x_min,x_min],[y_min,y_min,y_max,y_max,y_min];color=:black,linewidth=2)
    return nothing
end

function add_observation_box!(ax,center,half_width)
    isnothing(center) && return nothing
    x0,y0 = center
    hx,hy = half_width
    lines!(ax,[x0-hx,x0+hx,x0+hx,x0-hx,x0-hx],
           [y0-hy,y0-hy,y0+hy,y0+hy,y0-hy];color=:magenta,linewidth=2,linestyle=:dash)
    return nothing
end

function plot_forward_fields(stokes,grid,it,t,figdir;
                             x_min=-0.5,x_max=0.5,y_min=-0.05,y_max=0.05)
    (;xci,xvi) = grid
    ni = size(stokes.P)
    Vx_v = @zeros(ni .+ 1...)
    Vy_v = @zeros(ni .+ 1...)
    velocity2vertex!(Vx_v,Vy_v,@velocity(stokes)...)
    Vx,Vy = Array(Vx_v),Array(Vy_v)
    P = Array(stokes.P)
    ηc = Array(stokes.viscosity.η)
    τxx,τyy = Array(stokes.τ.xx),Array(stokes.τ.yy)
    τxy_c = Array(stokes.τ.xy_c)
    τzz = @. -τxx-τyy
    τII = @. sqrt(0.5*(τxx^2+τyy^2+τzz^2)+τxy_c^2)

    f1 = Figure(size=(1900,950))
    Label(f1[0,1:8],@sprintf("Forward solution — step %d, t = %.5f",it,t);fontsize=24)
    specs = (("Vx",xvi,Vx,:vik,true),("Vy",xvi,Vy,:vik,true),
             ("P",xci,P,:vik,true),("ηc",xci,ηc,:batlow,false),
             ("τxx",xci,τxx,:vik,true),("τyy",xci,τyy,:vik,true),
             ("τxy,c",xci,τxy_c,:vik,true),("τII",xci,τII,:batlow,false))
    for (k,(title,x,A,cmap,div)) in enumerate(specs)
        row,col = divrem(k-1,4)
        ax = panel!(f1,row+1,col+1,title,x,A;colormap=cmap,diverging=div)
        add_shear_zone!(ax,x_min,x_max,y_min,y_max)
    end
    save_figure(figdir,"forward"*step_tag(it),f1)
    return nothing
end

function sensitivity_fields(stokes_ad,gradients,grid)
    (;xci,xvi) = grid
    dJdG = dropdims(sum(Array(gradients.G.center),dims=1),dims=1)
    return ((title=L"\frac{dJ}{d\eta}",x=xvi,field=Array(stokes_ad.viscosity.ηv)),
            (title=L"\frac{dJ}{dG}",x=xci,field=dJdG),
            (title=L"\frac{dJ}{d\rho}",x=xci,field=Array(stokes_ad.ρ)))
end

# One sensitivity panel: axis at fig[row,col], colorbar at fig[row,col+1]
function sensitivity_panel!(fig,row,col,x,y,A;
                            title,k=col,mode=:contour,
                            xlabel="x [m]",ylabel="y [m]",cbar_label="",vmax=nothing,
                            tickformat="%.2e",
                            x_min,x_max,y_min,y_max,
                            observation_center,observation_half_width,
                            nlevels=20,colormap=:balance)
    if isnothing(vmax)
        vals = filter(isfinite,vec(A))
        vmax = isempty(vals) ? 0.0 : maximum(abs,vals)
    end
    if vmax == 0.0
        @warn "Sensitivity field is zero" field=k
        vmax = 1.0
    end
    log_at(2) && @printf("Sensitivity %d: min=% .6e max=% .6e vmax=% .6e\n",k,minimum(A),maximum(A),vmax)

    ax = Axis(fig[row,col];title=title,xlabel=xlabel,ylabel=ylabel,aspect=DataAspect(),
              titlesize=28,xticklabelsize=18,yticklabelsize=18,xgridvisible=false,ygridvisible=false)
    if mode == :contour
        levels = collect(range(-vmax,vmax;length=nlevels+1))
        plt = contourf!(ax,x,y,A;levels=levels,colormap=colormap,extendlow=:auto,extendhigh=:auto)
        minimum(A) < 0 < maximum(A) &&
            contour!(ax,x,y,A;levels=[0.0],color=(:black,0.35),linewidth=0.8)
    else
        plt = heatmap!(ax,x,y,A;colormap=colormap,colorrange=(-vmax,vmax))
    end

    ticks = [-vmax,-vmax/2,0.0,vmax/2,vmax]
    Colorbar(fig[row,col+1],plt;width=20,ticklabelsize=17,labelsize=18,label=cbar_label,
             ticks=(ticks,[Printf.format(Printf.Format(tickformat),v) for v in ticks]))

    add_shear_zone!(ax,x_min,x_max,y_min,y_max)
    add_observation_box!(ax,observation_center,observation_half_width)
    xlims!(ax,first(x),last(x))
    ylims!(ax,first(y),last(y))
    return ax
end

function plot_sensitivities(stokes_ad,gradients,grid,figdir;
                            mode=:contour,
                            x_min=-500.0,x_max=500.0,
                            y_min=-50.0,y_max=50.0,
                            observation_center=(0.0,0.0),
                            observation_half_width=(100.0,100.0),
                            nlevels=20,
                            colormap=:balance,
                            step=nothing)
    mode in (:heatmap,:contour) || error("mode must be :heatmap or :contour")
    fields = sensitivity_fields(stokes_ad,gradients,grid)
    fig = Figure(size=(1800,550),fontsize=22,figure_padding=(10,25,10,10))
    for (k,f) in enumerate(fields)
        sensitivity_panel!(fig,1,2k-1,collect(f.x[1]),collect(f.x[2]),Array(f.field);
                           title=f.title,k=k,mode=mode,ylabel=k==1 ? "y [m]" : "",
                           x_min,x_max,y_min,y_max,
                           observation_center,observation_half_width,
                           nlevels,colormap)
    end
    colgap!(fig.layout,20)
    name = mode==:heatmap ? "sensitivity_heatmap" : "sensitivity_contours"
    save_figure(figdir,name*step_tag(step),fig)
    return fig
end

# log10 axes: one tick per decade, or 1–3 / 1–2–5 steps for less than two decades.
# Use as Axis(...; yscale=log10, yticks=decade_ticks, ytickformat=decade_labels)
function decade_ticks(vmin,vmax)
    (isfinite(vmin) && isfinite(vmax) && 0 < vmin < vmax) || return Float64[]
    e0,e1 = floor(Int,log10(vmin)),ceil(Int,log10(vmax))
    for ms in ((1,),(1,3),(1,2,5))
        v = [m*10.0^e for e in e0:e1 for m in ms if vmin <= m*10.0^e <= vmax]
        length(v) >= 2 && return v
    end
    return [vmin,vmax]
end
function decade_labels(vs)
    return map(vs) do v
        e = floor(Int,log10(v)+1e-9); m = v/10.0^e
        isapprox(m,1;atol=1e-6) ? latexstring("10^{$e}") :
            latexstring(@sprintf("%g",round(m;sigdigits=2))*"\times 10^{$e}")
    end
end

# J is in m/s: raw ∂J/∂p in m/s per unit of p, p ∂J/∂p plotted in cm/yr
const SENSITIVITY_UNITS = (
    raw = (η=L"\mathrm{m\,s^{-1}\,/\,(Pa\,s)}",G=L"\mathrm{m\,s^{-1}\,/\,Pa}",
           ρ=L"\mathrm{m\,s^{-1}\,/\,(kg\,m^{-3})}"),
    logarithmic = (η=L"\mathrm{cm\,yr^{-1}}",G=L"\mathrm{cm\,yr^{-1}}",ρ=L"\mathrm{cm\,yr^{-1}}"),
)

# Inset in the lower part of `ax` (layout cell `pos`) showing the window
# zoom = ((x0,x1),(y0,y1)) with the same contour levels; the window is outlined on `ax`.
function add_zoom_inset!(pos,ax,x,y,A;zoom,vmax,nlevels=20,colormap=:balance,
                         x_min,x_max,y_min,y_max)
    (x0,x1),(y0,y1) = zoom
    lines!(ax,[x0,x1,x1,x0,x0],[y0,y0,y1,y1,y0];color=:black,linewidth=1.2)
    w  = 0.62
    h  = w*(y1-y0)/(x1-x0)*(last(x)-first(x))/(last(y)-first(y))   # keeps the data aspect
    ai = Axis(pos;width=Relative(w),height=Relative(h),halign=0.55,valign=0.17,
              backgroundcolor=:white,xticklabelsize=13,yticklabelsize=13,
              xticks=[x0,0.0,x1],yticks=[y0,0.0,y1],
              xgridvisible=false,ygridvisible=false)
    translate!(ai.blockscene,0,0,150)                    # draw on top of the main axis
    levels = collect(range(-vmax,vmax;length=nlevels+1))
    contourf!(ai,x,y,A;levels=levels,colormap=colormap,extendlow=:auto,extendhigh=:auto)
    add_shear_zone!(ai,x_min,x_max,y_min,y_max)
    limits!(ai,x0,x1,y0,y1)
    return ai
end

# Step, t in years and, given τM = η_weak/G, t/τM above a figure
function time_stamp!(pos,step,t,τM)
    parts = String[]
    isnothing(step) || push!(parts,"\\mathrm{time step\\ $step}")
    if !isnothing(t)
        ts = @sprintf("t = %.3g\\ \\mathrm{yr}",t/SEC_PER_YR)
        isnothing(τM) || (ts *= @sprintf(" = %.3g\\,\\tau_M\\quad (\\tau_M = \\eta_\\mathrm{weak}/G = %.3g\\ \\mathrm{yr})",
                                         t/τM,τM/SEC_PER_YR))
        push!(parts,ts)
    end
    isempty(parts) || Label(pos,latexstring(join(parts,",\\ "));fontsize=22,halign=:left,tellwidth=false)
    return nothing
end

# η, G, ρ: ∂J/∂p (row 1) and p ∂J/∂p (row 2, one shared colour limit)
function plot_sensitivities_combined(stokes_ad,gradients,grid,logarithmic,coords,figdir;
                                     x_min=-500.0,x_max=500.0,
                                     y_min=-50.0,y_max=50.0,
                                     observation_center=(0.0,0.0),
                                     observation_half_width=(100.0,100.0),
                                     nlevels=20,
                                     colormap=:balance,
                                     step=nothing,t=nothing,τM=nothing,
                                     zoom=((-100.0,100.0),(-50.0,50.0)))

    raw = sensitivity_fields(stokes_ad,gradients,grid)
    ks  = (:η,:G,:ρ)                                  # same order as sensitivity_fields
    log_titles = (η=L"\eta\,\frac{dJ}{d\eta}",G=L"G\,\frac{dJ}{dG}",ρ=L"\rho\,\frac{dJ}{d\rho}")
    common = (x_min=x_min,x_max=x_max,y_min=y_min,y_max=y_max,
              observation_center=observation_center,
              observation_half_width=observation_half_width,
              nlevels=nlevels,colormap=colormap)

    fig = Figure(size=(1800,1100),fontsize=22,figure_padding=(10,25,10,10))
    time_stamp!(fig[0,1:6],step,t,τM)

    logcm = NamedTuple{ks}(map(p->Array(logarithmic[p]) .* MS_TO_CMYR,ks))
    vmax_log = maximum(maximum(abs,filter(isfinite,vec(logcm[p]));init=0.0) for p in ks)
    log_at(2) && @printf("Logarithmic sensitivities: shared color limit ±%.6e cm/yr\n",vmax_log)

    for (k,p) in enumerate(ks)
        f = raw[k]
        sensitivity_panel!(fig,1,2k-1,collect(f.x[1]),collect(f.x[2]),Array(f.field);
                           title=f.title,k=k,xlabel="",ylabel=k==1 ? "y [m]" : "",
                           cbar_label=SENSITIVITY_UNITS.raw[p],common...)
        x,y = coords[p]
        ax = sensitivity_panel!(fig,2,2k-1,collect(x),collect(y),logcm[p];
                                title=log_titles[p],k=k,ylabel=k==1 ? "y [m]" : "",
                                cbar_label=SENSITIVITY_UNITS.logarithmic[p],vmax=vmax_log,
                                tickformat="%.3g",common...)
        p === :η && !isnothing(zoom) &&
            add_zoom_inset!(fig[2,2k-1],ax,collect(x),collect(y),logcm[p];
                            zoom=zoom,vmax=vmax_log,nlevels=nlevels,colormap=colormap,
                            x_min=x_min,x_max=x_max,y_min=y_min,y_max=y_max)
    end

    colgap!(fig.layout,20)
    rowgap!(fig.layout,20)
    save_figure(figdir,"sensitivity_combined"*step_tag(step),fig)
    return fig
end

## Gradient descent diagnostics ---------------------------------------------------------
"""
    plot_gradient_descent_G(result; G_true=nothing, figdir=nothing, name="gradient_descent_G")

Misfit-over-G plot of `gradient_descent_G`, built from `result.history`.
- Solid line + filled markers: the accepted path (start point and every accepted step).
- Faint dashed segments + hollow markers: rejected backtracking trials, each drawn from the
  accepted point it was tried from.
`G_true` (optional) adds a vertical line at the synthetic-truth G. `frozen` (optional, from
`frozen_history_curves`) adds the misfit with τ_o frozen around every gradient point: only the
last time step sees G, which is the function the per-step adjoint differentiates. If `figdir`
is given the figure is saved there.
"""
function plot_gradient_descent_G(result; G_true=nothing, figdir=nothing, name="gradient_descent_G", yscale=identity,
                                 frozen=nothing)
    h = result.history
    # the start point (gd_iter = 0) carries accepted = false but belongs to the accepted path
    is_acc = [e.accepted for e in h]
    path   = h[is_acc]
    trials = h[.!is_acc]

    # base point of every rejected trial: the last accepted point before it
    base = Vector{eltype(h)}(undef, length(trials))
    last_acc = first(h); k = 0
    for (e, a) in zip(h, is_acc)
        a ? (last_acc = e) : (k += 1; base[k] = last_acc)
    end

    fig = Figure(size=(800,650),fontsize=20)
    ax = Axis(fig[1,1];xlabel="G [GPa]",ylabel="misfit",yscale=yscale,
              title="Gradient descent on G ($(length(path)-1) accepted, $(length(trials)) rejected)")

    c_acc, c_rej = :dodgerblue3, :darkorange2

    for (b, e) in zip(base, trials)       # rejected trial segments, from where they were tried
        lines!(ax,[b.G,e.G] ./ 1e9,[b.m,e.m];color=(c_rej,0.35),linestyle=:dash,linewidth=1.5)
    end
    scatter!(ax,[e.G for e in trials] ./ 1e9,[e.m for e in trials];color=(c_rej,0.15),
             strokecolor=(c_rej,0.6),strokewidth=1.5,markersize=12,
             label="rejected trial (m_try ≥ current m)")

    lines!(ax,[e.G for e in path] ./ 1e9,[e.m for e in path];color=c_acc,linewidth=2.5)
    scatter!(ax,[e.G for e in path] ./ 1e9,[e.m for e in path];color=c_acc,markersize=12,
             label="accepted")
    # misfit with τ_o frozen at each gradient point (drawn under the tangents)
    if !isnothing(frozen)
        for (k, c) in enumerate(frozen)
            lines!(ax,c.Gs ./ 1e9,c.ms;color=(:mediumpurple3,0.9),linewidth=2,
                   label = k == 1 ? "misfit, τ_o frozen (last step only)" : nothing)
        end
    end
    # local misfit slope at every accepted point: a short black tangent segment. The gradient at a
    # point is computed when the next iteration starts from it, so it sits in the history entries
    # of iteration gd_iter+1 (the final point has none).
    xs, ys = [e.G for e in path] ./ 1e9, [e.m for e in path]
    xspan = max(maximum(xs) - minimum(xs), eps())
    yspan = max(maximum(ys) - minimum(ys), eps())
    first_grad = false
    for e in path
        j = findfirst(t -> t.gd_iter == e.gd_iter + 1, h)
        isnothing(j) && continue
        slope = h[j].dmdG * 1e9                       # dm/dG per GPa, matching the x axis
        # half-width: 4 % of the G range, shortened so the segment never exceeds 8 % of the misfit
        # range (or 30 % of the point's own misfit, which keeps it positive on a log axis)
        dx = min(0.04 * xspan, slope == 0 ? Inf : min(0.08 * yspan, 0.3 * e.m) / abs(slope))
        lines!(ax,[e.G/1e9 - dx, e.G/1e9 + dx],[e.m - slope*dx, e.m + slope*dx];color=:black,linewidth=3,
               label = first_grad ? nothing : "gradient dm/dG")
        first_grad = true
    end

    scatter!(ax,[first(path).G] ./ 1e9,[first(path).m];color=:black,marker=:diamond,markersize=16,
             label="start (G₀)")
    scatter!(ax,[last(path).G] ./ 1e9,[last(path).m];color=:red,marker=:star5,markersize=20,
             label="final ($(get(result,:status,"")))")

    isnothing(G_true) || vlines!(ax,G_true/1e9;color=:gray40,linestyle=:dot,linewidth=2,label="G true")

    # compact legend under the axis: small font, 3 columns, row height sized to its content so the
    # main axis keeps (nearly) all of the figure height
    Legend(fig[2,1],ax;framevisible=false,labelsize=20,nbanks=3,orientation=:vertical,
           markersize=10,patchsize=(22,10),rowgap=0,colgap=18,padding=(0,0,0,0),
           tellheight=true,tellwidth=false)
    rowsize!(fig.layout,2,Auto(1))
    rowgap!(fig.layout,6)
    isnothing(figdir) || save_figure(figdir,"$(name)_$(nameof(yscale))",fig)   # e.g. gradient_descent_G_log.png
    return fig
end

## Viscosity (ηv) inversion diagnostics ------------------------------------------------
# All take the named tuple returned by `gradient_descent_η` and the model `grid`
# (`viscosity_grid(setup)` builds one without a model). ηv lives on the vertices: grid.xvi.
# Column nx+1 of ηv is the periodic seam (copy of column 1) and is left out of x-averages.

viscosity_grid(setup) = Geometry((setup.domain.nx,setup.domain.ny),(setup.domain.lx,setup.domain.ly);
                                 origin=setup.domain.origin)

# the start point and every accepted step / the rejected backtracking trials
ηpath(h)   = [e for e in h if e.accepted || e.gd_iter == 0]
ηtrials(h) = [e for e in h if !(e.accepted || e.gd_iter == 0)]

# x-average of log10 ηv (geometric mean), seam column excluded, plus the min/max over x
function ηprofile(ηv)
    L = log10.(ηv[1:end-1,:])
    return (; mean=vec(sum(L;dims=1)) ./ size(L,1),lo=vec(minimum(L;dims=1)),hi=vec(maximum(L;dims=1)))
end

# gradient of iteration i (taken at the accepted point of iteration i-1), as d ln m / d ln ηv
function ηgradient(h,i)
    j = findfirst(e -> e.gd_iter == i,h)
    isnothing(j) && return nothing
    m_base = ηpath(h)[findlast(e -> e.gd_iter < i,ηpath(h))].m
    return h[j].dmdlnηv ./ m_base
end

# iteration axis: integer ticks only (steps of 1, 2, 5, 10, ...)
const ITER_TICKS = WilkinsonTicks(6;k_min=3,Q=[(1.0,1.0),(5.0,0.9),(2.0,0.7)])

# log axis with the decade ticks defined above
logaxis(pos;kwargs...) = Axis(pos;yscale=log10,yticks=decade_ticks,ytickformat=decade_labels,kwargs...)

"""
    plot_gradient_descent_η(result; figdir=nothing, name="gradient_descent_eta")

Convergence of `gradient_descent_η`, against the gradient-descent iteration:
(a) misfit m/m₀, accepted path and rejected backtracking trials;
(b) largest accepted change max|Δ ln ηv| and the gradient max|d ln m/d ln ηv|;
(c) forward runs per iteration (backtracking trials), to judge the γ adaptation.
"""
function plot_gradient_descent_η(result;figdir=nothing,name="gradient_descent_eta")
    h = result.history
    path,trials = ηpath(h),ηtrials(h)
    m0 = first(h).m
    niter = maximum(e.gd_iter for e in h)
    c_acc,c_rej = :dodgerblue3,:darkorange2

    fig = Figure(size=(1700,560),fontsize=20)
    Label(fig[0,1:3],"Viscosity inversion: $(length(path)-1) accepted steps, $(length(trials)) rejected trials, " *
          "status $(get(result,:status,""))";fontsize=22)

    ax1 = logaxis(fig[1,1];xlabel="iteration",xticks=ITER_TICKS,ylabel=L"m/m_0",title="misfit")
    isempty(trials) || scatter!(ax1,[e.gd_iter for e in trials],[e.m/m0 for e in trials];
                                color=(c_rej,0.15),strokecolor=(c_rej,0.7),strokewidth=1.5,markersize=10,
                                label="rejected trial")
    scatterlines!(ax1,[e.gd_iter for e in path],[e.m/m0 for e in path];color=c_acc,markersize=10,
                  linewidth=2.5,label="accepted")
    axislegend(ax1;position=:rt,framevisible=false)

    ax2 = logaxis(fig[1,2];xlabel="iteration",xticks=ITER_TICKS,title="step size and gradient")
    its  = [e.gd_iter for e in path[2:end]]
    dlnη = [maximum(abs,log.(path[k+1].ηv ./ path[k].ηv)) for k in 1:length(path)-1]
    grad = [(i,maximum(abs,g)) for i in 1:niter for g in (ηgradient(h,i),) if !isnothing(g)]
    isempty(dlnη) || scatterlines!(ax2,its,dlnη;color=c_acc,markersize=10,linewidth=2,
                                   label=L"\max\,|\Delta \ln \eta_v|\ \mathrm{(accepted)}")
    scatterlines!(ax2,first.(grad),last.(grad);color=:black,marker=:utriangle,markersize=10,linewidth=2,
                  label=L"\max\,|\partial \ln m/\partial \ln \eta_v|")
    axislegend(ax2;position=:rt,framevisible=false)

    ax3 = Axis(fig[1,3];xlabel="iteration",ylabel="forward runs",title="line-search trials per iteration",
               xticks=ITER_TICKS,yticks=WilkinsonTicks(5;k_min=3))
    counts = [count(e -> e.gd_iter == i,h) for i in 1:niter]
    ok     = [any(e -> e.gd_iter == i && e.accepted,h) for i in 1:niter]
    barplot!(ax3,1:niter,counts;color=[a ? c_acc : :firebrick for a in ok],strokewidth=0)
    ylims!(ax3,0,nothing)

    colgap!(fig.layout,30)
    isnothing(figdir) || save_figure(figdir,name,fig)
    return fig
end

# map at fig[row,2col-1], colorbar at fig[row,2col] (`fig` may also be a GridLayout), shear zone and observation box drawn
function ηmap!(fig,row,col,xv,yv,A;title,colormap=:batlow,colorrange=crange(A),cbar_label="",
               x_min,x_max,y_min,y_max,observation_center,observation_half_width)
    ax = Axis(fig[row,2col-1];aspect=DataAspect(),title=title,xlabel="x [m]",ylabel=col==1 ? "y [m]" : "",
              titlesize=22,xgridvisible=false,ygridvisible=false)
    hm = heatmap!(ax,xv,yv,A;colormap=colormap,colorrange=colorrange)
    Colorbar(fig[row,2col],hm;width=18,ticklabelsize=15,label=cbar_label)
    add_shear_zone!(ax,x_min,x_max,y_min,y_max)
    add_observation_box!(ax,observation_center,observation_half_width)
    limits!(ax,first(xv),last(xv),first(yv),last(yv))
    col > 1 && hideydecorations!(ax;ticks=false,ticklabels=false)
    return ax
end

"""
    plot_viscosity_fields(result, grid; ηv_true=nothing, figdir=nothing, zone/box keywords...)

Maps of the ηv inversion. Row 1: log10 ηv of the start, the result and (given `ηv_true`) the
truth, on one colour scale. Row 2: log10(ηv/ηv_true) (the error; without `ηv_true` the change
log10(ηv/ηv₀)), the change from the start, and the last gradient d ln m/d ln ηv.
"""
function plot_viscosity_fields(result,grid;ηv_true=nothing,figdir=nothing,name="viscosity_inversion_fields",
                               x_min=-500.0,x_max=500.0,y_min=-50.0,y_max=50.0,
                               observation_center=(0.0,0.0),observation_half_width=(100.0,100.0))
    h = result.history
    xv,yv = collect.(grid.xvi)
    L0,L1 = log10.(first(h).ηv),log10.(result.ηv)
    Lt = isnothing(ηv_true) ? nothing : log10.(ηv_true)
    lims = crange(vcat(vec(L0),vec(L1),isnothing(Lt) ? Float64[] : vec(Lt)))
    kw = (; x_min,x_max,y_min,y_max,observation_center,observation_half_width)
    unit = L"\log_{10}\,\eta_v\ \mathrm{[Pa\,s]}"

    fig = Figure(size=(1800,1150),fontsize=20)
    Label(fig[0,1:6],@sprintf("Viscosity inversion — m/m₀ = %.3e after %d iterations (%s)",
                              result.misfit/first(h).m,maximum(e.gd_iter for e in h),get(result,:status,""));
          fontsize=24)
    ηmap!(fig,1,1,xv,yv,L0;title=L"\mathrm{start}\ \eta_v^{0}",colorrange=lims,cbar_label=unit,kw...)
    ηmap!(fig,1,2,xv,yv,L1;title=L"\mathrm{result}\ \eta_v",colorrange=lims,cbar_label=unit,kw...)
    isnothing(Lt) || ηmap!(fig,1,3,xv,yv,Lt;title=L"\mathrm{true}\ \eta_v^{\mathrm{true}}",colorrange=lims,
                           cbar_label=unit,kw...)

    if !isnothing(Lt)
        E = L1 .- Lt
        ηmap!(fig,2,1,xv,yv,E;title=L"\mathrm{error}\ \log_{10}(\eta_v/\eta_v^{\mathrm{true}})",
              colormap=:vik,colorrange=crange(E;diverging=true),kw...)
    end
    D = L1 .- L0
    ηmap!(fig,2,2,xv,yv,D;title=L"\mathrm{change}\ \log_{10}(\eta_v/\eta_v^{0})",
          colormap=:vik,colorrange=crange(D;diverging=true),kw...)
    i_last = maximum(e.gd_iter for e in h)
    g = i_last > 0 ? ηgradient(h,i_last) : nothing
    isnothing(g) || ηmap!(fig,2,3,xv,yv,g;title=latexstring("\\partial \\ln m/\\partial \\ln \\eta_v\\ \\mathrm{(iteration\\ $i_last)}"),
                          colormap=:balance,colorrange=crange(g;diverging=true),kw...)

    colgap!(fig.layout,15)
    isnothing(figdir) || save_figure(figdir,name,fig)
    return fig
end

# shaded weak zone and dashed observation-box limits on a profile axis (y vertical)
function add_profile_bands!(ax,y_min,y_max,observation_center,observation_half_width)
    hspan!(ax,y_min,y_max;color=(:gray,0.15))
    isnothing(observation_center) ||
        hlines!(ax,observation_center[2] .+ [-1,1] .* observation_half_width[2];color=:magenta,linestyle=:dash,linewidth=1.5)
    return nothing
end

"""
    plot_viscosity_profiles(result, grid; ηv_true=nothing, figdir=nothing, ...)

x-averaged log10 ηv against y for every accepted iterate (coloured by iteration), with the
start, the truth and the result's min/max over x. Grey band: weak zone; dashed: observation box.
"""
function plot_viscosity_profiles(result,grid;ηv_true=nothing,figdir=nothing,name="viscosity_inversion_profiles",
                                 y_min=-50.0,y_max=50.0,
                                 observation_center=(0.0,0.0),observation_half_width=(100.0,100.0),
                                 colormap=:viridis)
    path = ηpath(result.history)
    yv = collect(grid.xvi[2])
    its = [e.gd_iter for e in path]
    crng = (0,max(maximum(its),1))

    fig = Figure(size=(800,900),fontsize=20)
    ax = Axis(fig[1,1];xlabel=L"\log_{10}\,\eta_v\ \mathrm{[Pa\,s]}\ \ (x\mathrm{-average})",ylabel="y [m]",
              title="ηv profiles during the inversion")
    add_profile_bands!(ax,y_min,y_max,observation_center,observation_half_width)
    for e in path[2:end-1]
        lines!(ax,ηprofile(e.ηv).mean,yv;color=e.gd_iter,colormap=colormap,colorrange=crng,linewidth=1.2)
    end
    p0,p1 = ηprofile(first(path).ηv),ηprofile(result.ηv)
    lines!(ax,p0.mean,yv;color=:black,linestyle=:dash,linewidth=2.5,label="start")
    isnothing(ηv_true) || lines!(ax,ηprofile(ηv_true).mean,yv;color=(:gray30,0.8),linewidth=6,label="true")
    lines!(ax,p1.mean,yv;color=:red,linewidth=3,label="result")
    lines!(ax,p1.lo,yv;color=:red,linestyle=:dot,linewidth=1.5,label="result min/max over x")
    lines!(ax,p1.hi,yv;color=:red,linestyle=:dot,linewidth=1.5)
    ylims!(ax,first(yv),last(yv))
    Colorbar(fig[1,2];colormap=colormap,limits=crng,label="iteration")
    axislegend(ax;position=:rb,framevisible=false,labelsize=17)
    isnothing(figdir) || save_figure(figdir,name,fig)
    return fig
end

"""
    plot_viscosity_frames(result, grid; ηv_true=nothing, figdir, every=1, ...)

One frame per `every`-th accepted iterate (and the last): log10 ηv map next to its x-averaged
profile, on fixed axes, for an animation. Saved as viscosity_inversion_stepNNNN.png in
`figdir/viscosity_inversion/` (NNNN = gradient-descent iteration).
"""
function plot_viscosity_frames(result,grid;ηv_true=nothing,figdir,every=1,
                               x_min=-500.0,x_max=500.0,y_min=-50.0,y_max=50.0,
                               observation_center=(0.0,0.0),observation_half_width=(100.0,100.0))
    path = ηpath(result.history)
    xv,yv = collect.(grid.xvi)
    m0 = first(path).m
    Ls = [log10.(e.ηv) for e in path]
    isnothing(ηv_true) || push!(Ls,log10.(ηv_true))
    lims = crange(reduce(vcat,vec.(Ls)))
    kw = (; x_min,x_max,y_min,y_max,observation_center,observation_half_width)
    p0 = ηprofile(first(path).ηv).mean
    pt = isnothing(ηv_true) ? nothing : ηprofile(ηv_true).mean

    for (k,e) in enumerate(path)
        (k == length(path) || (k-1) % every == 0) || continue
        fig = Figure(size=(1400,700),fontsize=20)
        Label(fig[0,1:2],@sprintf("iteration %d — m/m₀ = %.3e",e.gd_iter,e.m/m0);fontsize=24)
        map_layout = fig[1,1] = GridLayout()
        ηmap!(map_layout,1,1,xv,yv,log10.(e.ηv);title=L"\log_{10}\,\eta_v",colorrange=lims,
              cbar_label=L"\mathrm{[Pa\,s]}",kw...)
        ax = Axis(fig[1,2];xlabel=L"\log_{10}\,\eta_v\ (x\mathrm{-average})",ylabel="",title="profile")
        add_profile_bands!(ax,y_min,y_max,observation_center,observation_half_width)
        lines!(ax,p0,yv;color=:black,linestyle=:dash,linewidth=2,label="start")
        isnothing(pt) || lines!(ax,pt,yv;color=(:gray30,0.8),linewidth=6,label="true")
        lines!(ax,ηprofile(e.ηv).mean,yv;color=:red,linewidth=3,label="current")
        limits!(ax,lims...,first(yv),last(yv))
        axislegend(ax;position=:rb,framevisible=false,labelsize=16)
        colsize!(fig.layout,2,Relative(0.33))
        save_figure(figdir,"viscosity_inversion"*step_tag(e.gd_iter),fig)
    end
    return nothing
end

"""
    plot_Vx_fit(result, grid, Vx_obs; figdir=nothing, ...)

Data fit of the ηv inversion, from `result.Vx0`, `result.Vx` and `result.mask`.
(a) x-averaged Vx against y: data, start and result, observed y range dashed;
(b, c) Vx - Vx_obs of the start and the result, observation box dashed. Vx in cm/yr.
"""
function plot_Vx_fit(result,grid,Vx_obs;figdir=nothing,name="viscosity_inversion_Vx_fit",
                     x_min=-500.0,x_max=500.0,y_min=-50.0,y_max=50.0,
                     observation_center=(0.0,0.0),observation_half_width=(100.0,100.0))
    xs,ys = collect.(grid.xi_vel[1])
    cm(V) = Array(V) .* MS_TO_CMYR
    Vobs,V0,V1 = cm(Vx_obs),cm(result.Vx0),cm(result.Vx)
    xavg(V) = vec(sum(V[1:end-1,:];dims=1)) ./ (size(V,1)-1)     # seam column excluded
    kw = (; x_min,x_max,y_min,y_max,observation_center,observation_half_width)

    fig = Figure(size=(1800,650),fontsize=20)
    ax = Axis(fig[1,1];xlabel=L"V_x\ \mathrm{[cm/yr]}\ \ (x\mathrm{-average})",ylabel="y [m]",title="velocity profile")
    add_profile_bands!(ax,y_min,y_max,observation_center,observation_half_width)
    lines!(ax,xavg(Vobs),ys;color=(:gray30,0.8),linewidth=6,label="data")
    lines!(ax,xavg(V0),ys;color=:black,linestyle=:dash,linewidth=2,label="start")
    lines!(ax,xavg(V1),ys;color=:red,linewidth=2.5,label="result")
    ylims!(ax,first(ys),last(ys))
    axislegend(ax;position=:rb,framevisible=false,labelsize=17)

    R0,R1 = V0 .- Vobs,V1 .- Vobs
    lab = L"\mathrm{[cm/yr]}"
    maps = fig[1,2] = GridLayout()
    ηmap!(maps,1,1,xs,ys,R0;title=L"\mathrm{start}\ V_x - V_x^{\mathrm{obs}}",colormap=:vik,
          colorrange=crange(R0;diverging=true),cbar_label=lab,kw...)
    ηmap!(maps,1,2,xs,ys,R1;title=L"\mathrm{result}\ V_x - V_x^{\mathrm{obs}}",colormap=:vik,
          colorrange=crange(R1;diverging=true),cbar_label=lab,kw...)
    colsize!(fig.layout,1,Relative(0.22))
    isnothing(figdir) || save_figure(figdir,name,fig)
    return fig
end

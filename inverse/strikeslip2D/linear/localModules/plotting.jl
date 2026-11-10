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
    "sensitivity_relative_shared"=>"relative_shared","sensitivity_max_relative"=>"max_relative")
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

### A Pluto.jl notebook ###
# v0.20.10

using Markdown
using InteractiveUtils

# ╔═╡ b3bc56eb-43a7-4736-bd66-704529911d60
begin
	using Pkg; Pkg.activate()
	using DrWatson
	md"Using DrWatson to ensure reproducibility between different machines ..."
end

# ╔═╡ 9071763e-f6ad-4468-ae48-369307a85263
begin
	@quickactivate "2024GLConvIso"
	using Dates
	using DelimitedFiles
	using ERA5Reanalysis
	using NASAPrecipitation
	using NCDatasets
	using Printf
	using Statistics
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("ultraplot")

	include(srcdir("common.jl"))
	
	md"Loading modules for the 2024GLConvIso project..."
end

# ╔═╡ a8431d46-46fe-11ec-2b8d-e39caffdabec
md"
# Figure 3b. Seasonal Variability of Rainfall over Observational Stations

In this notebook, we explicitly compare the rainfall measured by the 3 stations we have, to GPM Rainfall data at the nearest gridpoint corresponding to these stations.
"

# ╔═╡ a27a0d50-d460-4908-b0a6-c4bdf6b19d21
function calculatebufferweights(shiftsteps)

    buffer = Int(ceil((shiftsteps-1)/2))
    weights = ones(buffer*2+1)
    if buffer >= (shiftsteps/2)
        weights[1] = 0.5
        weights[end] = 0.5
    end
    weights /= shiftsteps
    return buffer,weights

end

# ╔═╡ 64c5c24e-6c15-444f-abc5-0e04afe562e8
function smooth(data::AbstractVector,days)

        buffer,weights = calculatebufferweights(days)

        ndt = length(data)
        ndata = fill(NaN,ndt)
        smth  = zeros(1+buffer*2)

        for ii in (1+buffer) : (ndt-buffer)

                for ismth = 0 : (buffer*2)
                        smth[ismth+1] = data[ii+ismth-buffer] * weights[ismth+1]
                end
                ndata[ii] = sum(smth)

        end

        return ndata

end

# ╔═╡ 507ed308-9a21-4f74-92bb-1020acc7dacd
begin
	infody  = stninfody(); ndystn = size(infody,1)
	md"Loading station location information ..."
end

# ╔═╡ b74efe40-6288-4b27-b757-5d0771f2552e
begin
	pplt.close()
	fig,axs = pplt.subplots(ncols=3,nrows=4,aspect=3,axwidth=1.5,hspace=1,wspace=1)

	for istn = 1 : 12

		stnstr = @sprintf("%02d",istn)
		geo = GeoRegion("OTREC_wrf_stn$stnstr",path=srcdir())

		wrfds = NCDataset(datadir("wrf","processed","$(geo.ID)-rain-daily-20190801_20201231-smooth_30days.nc"))
		wprcp = wrfds["RAINNC"][:]
		close(wrfds)

		gpmds = NCDataset(datadir("wrf","processed","$(geo.ID)-gpmrain-20190801_20201231.nc"))
		gprcp = dropdims(mean(reshape(gpmds["precipitation"][:],24,:),dims=1),dims=1) * 86400
		gprcp = smooth(gprcp,30)
		close(gpmds)

		erads = NCDataset(datadir("wrf","processed","$(geo.ID)-erarain-20190801_20201231.nc"))
		eprcp = dropdims(sum(reshape(erads["precipitation"][:],24,:),dims=1),dims=1)
		eprcp = smooth(eprcp,30)
		close(erads)

		if istn != 6
			axs[istn].plot(Date(2019,8,1):Day(1):Date(2020,12,31),wprcp)
			axs[istn].plot(Date(2019,8,1):Day(1):Date(2020,12,31),gprcp)
			axs[istn].plot(Date(2019,8,1):Day(1):Date(2020,12,31),eprcp)
			axs[istn].format(ultitle="$(infody[istn,1])")
		else
			axs[istn].plot(Date(2019,8,1):Day(1):Date(2020,12,31),wprcp,legend="r",legend_kw=Dict("frame"=>false,"ncol"=>1),label="WRF")
			axs[istn].plot(Date(2019,8,1):Day(1):Date(2020,12,31),gprcp,legend="r",legend_kw=Dict("frame"=>false,"ncol"=>1),label="IMERGv7")
			axs[istn].plot(Date(2019,8,1):Day(1):Date(2020,12,31),eprcp,legend="r",legend_kw=Dict("frame"=>false,"ncol"=>1),label="ERA5")
			axs[istn].format(ultitle="$(infody[istn,1])")
		end

	end
	

	for ax in axs
		ax.format(
			xlabel = "Date",
			ylabel = L"Precipitation / mm day$^{-1}$", ylim = (0,60),
			suptitle = "(b) 30-Day Smoothed Rainfall Rate"
		)
	end
	
	fig.savefig(projectdir("figures","fig3b-seasonal.png"),transparent=false,dpi=400)
	load(projectdir("figures","fig3b-seasonal.png"))
end

# ╔═╡ Cell order:
# ╠═a8431d46-46fe-11ec-2b8d-e39caffdabec
# ╟─b3bc56eb-43a7-4736-bd66-704529911d60
# ╟─9071763e-f6ad-4468-ae48-369307a85263
# ╟─a27a0d50-d460-4908-b0a6-c4bdf6b19d21
# ╟─64c5c24e-6c15-444f-abc5-0e04afe562e8
# ╟─507ed308-9a21-4f74-92bb-1020acc7dacd
# ╟─b74efe40-6288-4b27-b757-5d0771f2552e

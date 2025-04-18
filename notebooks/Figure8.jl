### A Pluto.jl notebook ###
# v0.20.5

using Markdown
using InteractiveUtils

# ╔═╡ 1d7446ba-f464-44df-89e2-ae2a5726e849
begin
	using Pkg; Pkg.activate()
	using DrWatson
	md"Using DrWatson in order to ensure reproducibility between different machines ..."
end

# ╔═╡ ab294fae-101f-4587-a2f4-7d72254dd421
begin
	@quickactivate "2024GLConvIso"
	using DelimitedFiles
	using GeoRegions
	using NCDatasets
	using PlutoUI
	using Printf
	using StatsBase

	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("ultraplot")

	include(srcdir("common.jl"))

	md"Loading modules for the ConvectionIsotopes project..."
end

# ╔═╡ fa2f8740-f813-11ec-00e1-112e2dfacda7
md"
# 01d. Creating GeoRegions for Stations

In this notebook, we define additional GeoRegions of interest for plotting and for analysis based on WRF modelling output and as necessary for figures.
"

# ╔═╡ ccdf14a6-b200-42db-8455-0ea4eeb5ae2d
TableOfContents()

# ╔═╡ 582476a4-e280-458b-ac4c-7c681ff96a74
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

# ╔═╡ 6aff97ec-0bd3-4d84-9d5c-93393941ca4e
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

# ╔═╡ 10d1c691-00a7-47de-a8ca-8debcd3346c1
function extractbudget(geoname,days)

	dystr  = "daily-20190801_20201231-smooth_$(@sprintf("%02d",days))days"

	ds = NCDataset(datadir(
		"wrf","processed",
		"$geoname-QBUDGET-20190801_20201231.nc"
	))
	prcp = smooth(dropdims(sum(reshape(ds["P"][:],24,:),dims=1),dims=1),days)
	evap = smooth(dropdims(sum(reshape(ds["E"][:],24,:),dims=1),dims=1),days)
	close(ds)

	ds = NCDataset(datadir(
		"wrf","processed",
		"$geoname-∇decompose-20190801_20201231.nc"
	))
	advc = smooth(dropdims(mean(reshape(ds["ADV"][:],24,:),dims=1),dims=1),days) * 86400
	divg = smooth(dropdims(mean(reshape(ds["DIV"][:],24,:),dims=1),dims=1),days) * 86400
	close(ds)

	dsp = NCDataset(datadir(
		"wrf","processed",
		"$geoname-p_wwgt-$dystr.nc"
	))
	pwgt = dsp["p_wwgt"][:] / 100
	pwgt[(pwgt.>1000).|(pwgt.<0)] .= NaN
	pwgt = dsp["σ_wwgt"][:]
	pwgt[(pwgt.>1).|(pwgt.<0)] .= NaN
	close(dsp)

	return pwgt,prcp,evap,advc,divg
	
end

# ╔═╡ d8558ea0-a753-4693-8dbe-2dc9ea86b5a0
function extract(geoname,days)

	dystr  = "daily-20190801_20201231-smooth_$(@sprintf("%02d",days))days"

	ds = NCDataset(datadir("wrf","processed","$geoname-cdhodq-$(dystr).nc"))
	c0 = ds["cdO18dH2O"][1,:] .- 1
	c1 = ds["cdO18dH2O"][2,:] * 1e6
	close(ds)

	return c0,c1
	
end

# ╔═╡ 405dd308-1b04-4487-b1b2-86ff17459167
function plotdeplete(
	axes,ii;
	nID, days=0, prfx = "", cinfo = false
)

	IDplt  = 1 : nID
	preplt = 400:50:1000; prevec = (preplt[1:(end-1)] .+ preplt[2:end]) ./ 2
	npre = length(prevec)
	δmat = zeros(nID,npre)

	for stn in 1 : nID
		stnstr = @sprintf("%02d",stn)
		geoname = "OTREC_wrf_$(prfx)$stnstr"
		c0,c1 = extract(geoname,days)
		pwgt,prcp,evap,advc,divg = extractbudget(geoname,days)
		it = ((prcp.+advc.-evap).>2.5) .& (.!isnan.(pwgt))
		μc0 = mean(c0[it])
		μc1 = mean(c1[it])
		δmat[stn,:] = μc1 * 1e-6 * (prevec*1e2 .- 1e5) .+ μc0
	end

	lvls = -25:2.5:0
	c = axes[ii].pcolormesh(IDplt,prevec,δmat'.*1000,extend="both",levels=lvls,cmap="viridis")

	if cinfo
		return c
	else
		return
	end

end

# ╔═╡ deaed5af-5700-418f-a6f1-05e0c0637d75
begin
	pplt.close(); fig,axs = pplt.subplots(ncols=3,aspect=1.5,axwidth=1.5)

	cbar =
	plotdeplete(axs,1,nID=25,days=7,cinfo=true,prfx="ITCZ")
	plotdeplete(axs,2,nID=25,days=7,cinfo=true,prfx="CrossITCZ")
	plotdeplete(axs,3,nID=25,days=7,cinfo=true,prfx="PAC2ATL")

	axs[1].format(ylim=(1000,400),title="ITCZ",xlabel="Region Number",ylabel=L"$p_{q\omega}$",suptitle=L"Replicating 7-day $\delta^{18}$O using Coefficients from Linear Fit")
	axs[2].format(ylim=(1000,400),title="CrossITCZ",xlabel="Region Number")
	axs[3].format(ylim=(1000,400),title="PAC2ATL",xlabel="Region Number")

	fig.colorbar(cbar,label=L"$\delta^{18}$O",locator=-25:5:0)
	fig.savefig(projectdir("figures","fig8-idealizedreplication.png"),transparent=false,dpi=400)
	load(projectdir("figures","fig8-idealizedreplication.png"))
end

# ╔═╡ Cell order:
# ╟─fa2f8740-f813-11ec-00e1-112e2dfacda7
# ╟─1d7446ba-f464-44df-89e2-ae2a5726e849
# ╟─ab294fae-101f-4587-a2f4-7d72254dd421
# ╟─ccdf14a6-b200-42db-8455-0ea4eeb5ae2d
# ╟─582476a4-e280-458b-ac4c-7c681ff96a74
# ╟─6aff97ec-0bd3-4d84-9d5c-93393941ca4e
# ╟─10d1c691-00a7-47de-a8ca-8debcd3346c1
# ╠═d8558ea0-a753-4693-8dbe-2dc9ea86b5a0
# ╠═405dd308-1b04-4487-b1b2-86ff17459167
# ╠═deaed5af-5700-418f-a6f1-05e0c0637d75

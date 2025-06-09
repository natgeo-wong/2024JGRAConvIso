### A Pluto.jl notebook ###
# v0.20.10

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
	@quickactivate "2024JGRAConvIso"
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

	md"Loading modules for the 2024JGRAConvIso project..."
end

# ╔═╡ fa2f8740-f813-11ec-00e1-112e2dfacda7
md"
# 10. Idealized Replications of Depletion Values for a given $p_{q\omega,bl}$
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
		"$geoname-p_wwgt3-$dystr.nc"
	))
	pwgt = dsp["p_wwgt"][:] / 100
	pbl  = dsp["P"][end,:] / 100 * 0.8; pbl[pbl.>800] .= 800
	pwgt = pwgt .+ pbl
	pwgt[(pwgt.>1000).|(pwgt.<0)] .= NaN
	close(dsp)

	return pwgt,mean(pbl[.!isnan.(pbl)]),prcp,evap,advc,divg
	
end

# ╔═╡ d8558ea0-a753-4693-8dbe-2dc9ea86b5a0
function extract(geoname,days)

	dystr  = "daily-20190801_20201231-smooth_$(@sprintf("%02d",days))days"

	ds = NCDataset(datadir("wrf","processed","$geoname-cdhodq-$(dystr).nc"))
	c0 = ds["cdO18dH2O"][1,:] .- 1
	c1 = ds["cdO18dH2O"][2,:]
	close(ds)

	return c0*1e3,c1*1e6
	
end

# ╔═╡ 405dd308-1b04-4487-b1b2-86ff17459167
function plotdeplete(
	axes,ii;
	nID, days=0, prfx = "", cinfo = false, pmin, pmax
)

	IDplt  = 1 : nID
	preplt = 300:50:900; prevec = (preplt[1:(end-1)] .+ preplt[2:end]) ./ 2
	npre = length(prevec)
	δmat = zeros(nID,npre) * NaN

	for stn in 1 : nID
		stnstr = @sprintf("%02d",stn)
		geoname = "OTREC_wrf_$(prfx)$stnstr"
		c0,c1 = extract(geoname,days)
		pwgt,μpbl,prcp,evap,advc,divg = extractbudget(geoname,days)
		it = ((prcp.+advc.-evap).>pmin) .& ((prcp.+advc.-evap).<pmax) .& (.!isnan.(pwgt))
		if sum(it) > 10
			c0 = c0[it]; c0 = c0[(c0.>=-20).&(c0.<=0)]; μc0 = mean(c0)
			c1 = c1[it]; c1 = c1[(c1.>=0).&(c1.<=1)]; μc1 = mean(c1)
			δmat[stn,:] = μc1 * 1e-6 * (prevec*1e2 .- μpbl*100) .+ μc0/1e3
			δmat[stn,prevec.>=μpbl] .= μc0/1e3
		end
	end

	lvls = -25:-10
	c = axes[ii].pcolormesh(IDplt,prevec,δmat'.*1000,extend="both",levels=lvls,cmap="viridis")

	if cinfo
		return c
	else
		return
	end

end

# ╔═╡ deaed5af-5700-418f-a6f1-05e0c0637d75
begin
	pplt.close(); fig,axs = pplt.subplots(ncols=3,nrows=3,aspect=2,axwidth=1.5,hspace=1)

	cbar =
	plotdeplete(axs,1,nID=25,days=7,pmin=2.5,pmax=15,cinfo=true,prfx="ITCZ")
	plotdeplete(axs,2,nID=25,days=7,pmin=2.5,pmax=15,cinfo=true,prfx="CrossITCZ")
	plotdeplete(axs,3,nID=25,days=7,pmin=2.5,pmax=15,cinfo=true,prfx="PAC2ATL")
	plotdeplete(axs,4,nID=25,days=7,pmin=15,pmax=30,cinfo=true,prfx="ITCZ")
	plotdeplete(axs,5,nID=25,days=7,pmin=15,pmax=30,cinfo=true,prfx="CrossITCZ")
	plotdeplete(axs,6,nID=25,days=7,pmin=15,pmax=30,cinfo=true,prfx="PAC2ATL")
	plotdeplete(axs,7,nID=25,days=7,pmin=30,pmax=100,cinfo=true,prfx="ITCZ")
	plotdeplete(axs,8,nID=25,days=7,pmin=30,pmax=100,cinfo=true,prfx="CrossITCZ")
	plotdeplete(axs,9,nID=25,days=7,pmin=30,pmax=100,cinfo=true,prfx="PAC2ATL")

	axs[1].format(title="ITCZ",xlabel="Region Number",suptitle=L"Replicating 7-day $\delta^{18}$O using Coefficients from Linear Fit")
	axs[2].format(title="CrossITCZ",xlabel="Region Number")
	axs[3].format(title="PAC2ATL",xlabel="Region Number")

	for ax in axs
		ax.format(ylim=(900,350),ylocator=[500,750],ylabel=L"$p_{q\omega,bl}$ / hPa",leftlabels=["2.5 < P < 15","15 < P < 30","P > 30"])
	end

	fig.colorbar(cbar,label=L"$\delta^{18}$O",locator=-25:5:0)
	fig.savefig(projectdir("figures","fig10-idealizedreplication.png"),transparent=false,dpi=400)
	load(projectdir("figures","fig10-idealizedreplication.png"))
end

# ╔═╡ Cell order:
# ╟─fa2f8740-f813-11ec-00e1-112e2dfacda7
# ╟─1d7446ba-f464-44df-89e2-ae2a5726e849
# ╟─ab294fae-101f-4587-a2f4-7d72254dd421
# ╟─ccdf14a6-b200-42db-8455-0ea4eeb5ae2d
# ╟─582476a4-e280-458b-ac4c-7c681ff96a74
# ╟─6aff97ec-0bd3-4d84-9d5c-93393941ca4e
# ╟─10d1c691-00a7-47de-a8ca-8debcd3346c1
# ╟─d8558ea0-a753-4693-8dbe-2dc9ea86b5a0
# ╟─405dd308-1b04-4487-b1b2-86ff17459167
# ╟─deaed5af-5700-418f-a6f1-05e0c0637d75

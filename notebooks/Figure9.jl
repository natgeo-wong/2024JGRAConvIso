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
# 09. Binning the coefficients for $\partial_pq_h/\partial_pq$
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

	return pwgt,prcp,evap,advc,divg
	
end

# ╔═╡ d8558ea0-a753-4693-8dbe-2dc9ea86b5a0
function extract(geoname,days)

	dystr  = "daily-20190801_20201231-smooth_$(@sprintf("%02d",days))days"

	ds = NCDataset(datadir("wrf","processed","$geoname-cdhodq-$(dystr).nc"))
	c0 = ds["cdO18dH2O"][1,:] .- 1
	c1 = ds["cdO18dH2O"][2,:]
	close(ds)

	return c0*1000, c1*1e6
	
end

# ╔═╡ c793412d-71b6-4f2c-a9f4-15da6ec039e4
function plotcdqdp(
	axes,ii;
	nID, days=0, prfx = "", cinfo = false, pmin, pmax
)

	c0edge = -21 : 1
	c1edge = 0 : 0.05 : 1
	binc0 = zeros(length(c0edge)-1,nID)
	binc1 = zeros(length(c1edge)-1,nID)
	c0plt = (c0edge[1:(end-1)] .+ c0edge[2:end])/2
	c1plt = (c1edge[1:(end-1)] .+ c1edge[2:end])/2
	μc0 = zeros(nID)*NaN
	μc1 = zeros(nID)*NaN
	IDplt = 1 : nID

	for stn in 1 : nID
		stnstr = @sprintf("%02d",stn)
		geoname = "OTREC_wrf_$(prfx)$stnstr"
		c0,c1 = extract(geoname,days)
		pwgt,prcp,evap,advc,divg = extractbudget(geoname,days)
		it = ((prcp.+advc.-evap).>pmin) .& ((prcp.+advc.-evap).<pmax) .&  (.!isnan.(pwgt))
		binc0[:,stn] += fit(Histogram,c0[it],c0edge).weights
		binc1[:,stn] += fit(Histogram,c1[it],c1edge).weights
		c0 = c0[it]; c0 = c0[(c0.>=-20).&(c0.<=0)]
		c1 = c1[it]; c1 = c1[(c1.>=0).&(c1.<=1)]
		if sum(it) > 10
			μc0[stn] = mean(c0)
			μc1[stn] = mean(c1)
		end
	end

	lvls = 1:20
	c1 = 
	axes[ii+0].pcolormesh(IDplt,c0plt,binc0,extend="both",levels=lvls)
	axes[ii+9].pcolormesh(IDplt,c1plt,binc1,extend="both",levels=lvls)
	
	axes[ii+0].plot(IDplt,μc0)
	axes[ii+9].plot(IDplt,μc1)

	if cinfo
		return c1
	else
		return
	end

end

# ╔═╡ 30424aa0-cc38-4f50-8eb6-efd4f6c4c9c4
begin
	pplt.close(); fig,axs = pplt.subplots(nrows=6,ncols=3,aspect=1.5,axwidth=1.5,sharey=0,hspace=[1,1,2,1,1])

	c1 =
	plotcdqdp(axs,1,nID=25,pmin=2.5,pmax=15,days=7,cinfo=true,prfx="ITCZ")
	plotcdqdp(axs,2,nID=25,pmin=2.5,pmax=15,days=7,cinfo=true,prfx="CrossITCZ")
	plotcdqdp(axs,3,nID=25,pmin=2.5,pmax=15,days=7,cinfo=true,prfx="PAC2ATL")
	plotcdqdp(axs,4,nID=25,pmin=15,pmax=30, days=7,cinfo=true,prfx="ITCZ")
	plotcdqdp(axs,5,nID=25,pmin=15,pmax=30, days=7,cinfo=true,prfx="CrossITCZ")
	plotcdqdp(axs,6,nID=25,pmin=15,pmax=30, days=7,cinfo=true,prfx="PAC2ATL")
	plotcdqdp(axs,7,nID=25,pmin=30,pmax=100,days=7,cinfo=true,prfx="ITCZ")
	plotcdqdp(axs,8,nID=25,pmin=30,pmax=100,days=7,cinfo=true,prfx="CrossITCZ")
	plotcdqdp(axs,9,nID=25,pmin=30,pmax=100,days=7,cinfo=true,prfx="PAC2ATL")

	axs[1].format(title="ITCZ",xlabel="Region Number")
	axs[2].format(title="CrossITCZ",xlabel="Region Number")
	axs[3].format(title="PAC2ATL",xlabel="Region Number")
	axs[6].format(suptitle=L"Coefficients for Linear Fit on 7-day $\delta^{18}$O")

	axs[4].format(ylabel=L"$c_0 - 1$ / $\perthousand$")
	axs[13].format(ylabel=L"$c_1$ / 10$^{-1}$ $\perthousand$ hPa$^{-1}$")
	
	for ii = 1 : 9
		axs[ii].format(ylim=(-21,1),leftlabels=["2.5 < P < 15","15 < P < 30","P > 30","2.5 < P < 15","15 < P < 30","P > 30"])
		axs[ii+9].format(ylim=(0,1))
	end

	for ii = [2,3,5,6,8,9,11,12,14,15,17,18]
		axs[ii].format(yticklabels="none")
	end

	fig.colorbar(c1,label="Number of Observations",length=0.5)
	fig.savefig(projectdir("figures","fig9-bincdhdq.png"),transparent=false,dpi=400)
	load(projectdir("figures","fig9-bincdhdq.png"))
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
# ╟─c793412d-71b6-4f2c-a9f4-15da6ec039e4
# ╠═30424aa0-cc38-4f50-8eb6-efd4f6c4c9c4

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

# ╔═╡ fa426522-25e6-4b80-80f1-40d74d02e7fe
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

# ╔═╡ ec7a951c-bd32-40b5-9cc1-5326dec119bf
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
function extract(geoname,days)

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
	pwgt1 = dsp["p_wwgt"][:] / 100
	# pwgto[(pwgto.>1000).|(pwgto.<0)] .= NaN
	close(dsp)

	dsp = NCDataset(datadir(
		"wrf","processed",
		"$geoname-p_wwgt2-$dystr.nc"
	))
	pwgt2 = dsp["p_wwgt"][:] / 100
	# pwgto[(pwgto.>1000).|(pwgto.<0)] .= NaN
	close(dsp)

	dsp = NCDataset(datadir(
		"wrf","processed",
		"$geoname-p_wwgt3-$dystr.nc"
	))
	pwgt3 = dsp["p_wwgt"][:] / 100
	psfc  = dsp["P"][end,:] / 100; pbl = 0.8 * psfc; pbl[pbl.>800] .= 800
	# pwgtn[(pwgtn.>1000).|(pwgtn.<0)] .= NaN
	close(dsp)

	return pwgt1,pwgt2,pwgt3.+pbl,prcp,evap,advc,divg
	
end

# ╔═╡ 405dd308-1b04-4487-b1b2-86ff17459167
function plotnvo(
	axes; days=0
)

	pω = 500 : 10 : 1000
	pqω = 500 : 10 : 1000
	pqωbl = 500 : 10 : 1000

	binmat1 = zeros(length(pω) -1,length(pω) -1)
	binmat2 = zeros(length(pω) -1,length(pω) -1)
	binmat3 = zeros(length(pω) -1,length(pω) -1)
	
	for stn in 1 : 25
		stnstr = @sprintf("%02d",stn)
		geoname = "OTREC_wrf_ITCZ$stnstr"
		pwgt1,pwgt2,pwgt3,prcp,evap,advc,divg = extract(geoname,days)
		it = ((prcp.+advc.-evap).>2.5) .& (.!isnan.(pwgt1)).& (.!isnan.(pwgt2)) .& (.!isnan.(pwgt3))
		binmat1[:,:] += fit(Histogram,(pwgt2[it],pwgt1[it]),(pqω,pω)).weights
		binmat2[:,:] += fit(Histogram,(pwgt3[it],pwgt1[it]),(pqωbl,pω)).weights
		binmat3[:,:] += fit(Histogram,(pwgt2[it],pwgt3[it]),(pqω,pqωbl)).weights
		geoname = "OTREC_wrf_CrossITCZ$stnstr"
		pwgt1,pwgt2,pwgt3,prcp,evap,advc,divg = extract(geoname,days)
		it = ((prcp.+advc.-evap).>2.5) .& (.!isnan.(pwgt1)).& (.!isnan.(pwgt2)) .& (.!isnan.(pwgt3))
		binmat1[:,:] += fit(Histogram,(pwgt2[it],pwgt1[it]),(pqω,pω)).weights
		binmat2[:,:] += fit(Histogram,(pwgt3[it],pwgt1[it]),(pqωbl,pω)).weights
		binmat3[:,:] += fit(Histogram,(pwgt2[it],pwgt3[it]),(pqω,pqωbl)).weights
		geoname = "OTREC_wrf_PAC2ATL$stnstr"
		pwgt1,pwgt2,pwgt3,prcp,evap,advc,divg = extract(geoname,days)
		it = ((prcp.+advc.-evap).>2.5) .& (.!isnan.(pwgt1)).& (.!isnan.(pwgt2)) .& (.!isnan.(pwgt3))
		binmat1[:,:] += fit(Histogram,(pwgt2[it],pwgt1[it]),(pqω,pω)).weights
		binmat2[:,:] += fit(Histogram,(pwgt3[it],pwgt1[it]),(pqωbl,pω)).weights
		binmat3[:,:] += fit(Histogram,(pwgt2[it],pwgt3[it]),(pqω,pqωbl)).weights
	end

	lvls = [1,2,5,10,20,50,100,200,500,1000]
	c = axes[1].pcolormesh(levels=lvls,extend="both",pqω,pω,binmat1')
	axes[2].pcolormesh(levels=lvls,extend="both",pqωbl,pω,binmat2')
	axes[3].pcolormesh(levels=lvls,extend="both",pqω,pqωbl,binmat3')

	return c

end

# ╔═╡ deaed5af-5700-418f-a6f1-05e0c0637d75
begin
	pplt.close()
	fig,axs = pplt.subplots([[1,2],[3,0]],axwidth=1,sharex=0,sharey=0,wspace=1,hspace=1)

	c = plotnvo(axs,days=7)

	axs[1].format(xlabel=L"$p_{q\omega}$",ylabel=L"$p_\omega$",xtickloc="top")
	axs[2].format(xlabel=L"$p_{q\omega,bl}$",xtickloc="top",yticklabels="none")
	axs[3].format(ylabel=L"$p_{q\omega,bl}$")

	for ax in axs
		ax.plot([600,1000],[600,1000],c="grey3",lw=0.5)
		ax.format(
			xlim=(600,1000),xlocator=700:200:900,
			ylim=(600,1000),ylocator=700:200:900,
			suptitle=L"Comparing 7-day $p_\omega$, $p_{q\omega}$ and $p_{q\omega,bl}$"
		)
	end

	fig.colorbar(c,length=0.75)
	fig.savefig(projectdir("figures","figS8-compare.png"),transparent=false,dpi=400)
	load(projectdir("figures","figS8-compare.png"))
end

# ╔═╡ Cell order:
# ╟─fa2f8740-f813-11ec-00e1-112e2dfacda7
# ╟─1d7446ba-f464-44df-89e2-ae2a5726e849
# ╟─ab294fae-101f-4587-a2f4-7d72254dd421
# ╟─ccdf14a6-b200-42db-8455-0ea4eeb5ae2d
# ╟─fa426522-25e6-4b80-80f1-40d74d02e7fe
# ╟─ec7a951c-bd32-40b5-9cc1-5326dec119bf
# ╟─10d1c691-00a7-47de-a8ca-8debcd3346c1
# ╠═405dd308-1b04-4487-b1b2-86ff17459167
# ╟─deaed5af-5700-418f-a6f1-05e0c0637d75

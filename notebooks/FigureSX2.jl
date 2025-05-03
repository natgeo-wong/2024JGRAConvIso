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
	using ERA5Reanalysis
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

# ╔═╡ eea4a60c-274d-47e7-a655-6651ad867082
begin
	coast = readdlm(datadir("coast.cst"),comments=true)
	clon  = coast[:,1]
	clat  = coast[:,2]
	ii = (clon .<= -75) .& (clon .>= -90).& (clat .>= 0).& (clat .<= 15)
	clon[.!ii] .= NaN
	clat[.!ii] .= NaN
	md"Preloading coastline data"
end

# ╔═╡ b3df876e-bda6-4cc4-b305-f816168a4ae7
e5ds = ERA5Monthly(start=Date(2019),stop=Date(2020),path=datadir())

# ╔═╡ b21feaa1-89c2-4ef0-8744-8c8d6c95147e
evar = SingleVariable("sst",path=srcdir())

# ╔═╡ a4972591-26bd-4555-b96b-765f20051763
egeo = ERA5Region(GeoRegion("OTREC_wrf_d02",path=srcdir()))

# ╔═╡ 7e9229d8-08ff-4442-b028-684e52f4f668
elsd = getLandSea(e5ds,egeo); nlon = length(elsd.lon); nlat = length(elsd.lat);

# ╔═╡ b01fa482-2d15-470c-9e1d-35bf12fb7609
begin
	sst = zeros(nlon,nlat)
	ds = read(e5ds,evar,egeo,Date(2019))
	sst[:,:] += sum(ds["sst"][:,:,8:12],dims=3)
	close(ds)
	ds = read(e5ds,evar,egeo,Date(2020))
	sst[:,:] += sum(ds["sst"][:,:,:],dims=3)
	close(ds)
	sst /= 17
	md"Calculating climatology for SST"
end

# ╔═╡ acebf498-8cdf-475e-9609-f60e68d4a355
calcT2es(T::Real) = 611.657 * exp((2.5e6/461.5181) * (1/273.16 - 1/T))

# ╔═╡ 90522d32-25d7-45b7-aca5-a0c04da39432
calce2q(e::Real,p::Real) = e * 0.621981 / (p - 0.378019 * e)

# ╔═╡ 7c45e2c6-7d6a-42c2-afaf-4fd4d110b4bb
begin
	gds = NCDataset(datadir("wrf","grid.nc"))
	wln = gds["longitude"][:,:]
	wlt = gds["latitude"][:,:]
	wpb = gds["pressure_base"][:,:,1]
	close(gds)
end

# ╔═╡ 943902fb-6853-46d7-8a0c-c4ef3bcd8229
begin
	qds = NCDataset(datadir("wrf","3D","QVAPOR-daily-20190801_20201231.nc"))
	qsfc = dropdims(mean(qds["QVAPOR"][:,:,1,1:(end-1)],dims=3),dims=3) * 1000
	close(qds)
end

# ╔═╡ b3475222-d143-481f-85a5-82ee6272a6fa
begin
	pds = NCDataset(datadir("wrf","3D","P-daily-20190801_20201231.nc"))
	psfc = dropdims(mean(pds["P"][:,:,1,1:(end-1)],dims=3),dims=3) .+ wpb
	close(pds)
end

# ╔═╡ 29603238-6c17-4668-8594-d240a76a23b9
begin
	tds = NCDataset(datadir("wrf","3D","T-daily-20190801_20201231.nc"))
	tsfc = dropdims(mean(tds["T"][:,:,1,1:(end-1)],dims=3),dims=3) .+ 300
	tsfc = tsfc .* (100000 ./ psfc) .^ (287/1004)
	close(tds)
	esat = calcT2es.(tsfc)
	qsat = calce2q.(esat,psfc) * 1000
end

# ╔═╡ deaed5af-5700-418f-a6f1-05e0c0637d75
begin
	pplt.close(); fig,axs = pplt.subplots(ncols=2,axwidth=1.5)

	c1 = axs[1].pcolormesh(elsd.lon,elsd.lat,sst',extend="both",levels=297:0.5:302)
	c1 = axs[1].pcolormesh(wln,wlt,tsfc,extend="both",levels=297:0.5:302)
	c2 = axs[2].pcolormesh(wln,wlt,qsfc./qsat*100,extend="both",cmap="drywet",levels=75:95)

	axs[1].format(title=L"$\mu$(SST) / K")
	axs[2].format(title=L"$\mu(rh_{sfc})$ / %")

	axs[1].colorbar(c1)
	axs[2].colorbar(c2,locator=75:5:95)

	for ax in axs
		ax.plot(clon,clat,c="k",lw=0.5)
		ax.format(
			suptitle="Mean Climatology (Aug 2019 - Dec 2020)",
			xlim=(-90,-75),xlabel=L"Longitude / $\degree$",
			ylim=(0,15),ylabel=L"Latitude / $\degree$"
		)
	end

	fig.savefig(projectdir("figures","figSX2.png"),transparent=false,dpi=400)
	load(projectdir("figures","figSX2.png"))
end

# ╔═╡ Cell order:
# ╟─fa2f8740-f813-11ec-00e1-112e2dfacda7
# ╟─1d7446ba-f464-44df-89e2-ae2a5726e849
# ╟─ab294fae-101f-4587-a2f4-7d72254dd421
# ╟─ccdf14a6-b200-42db-8455-0ea4eeb5ae2d
# ╠═eea4a60c-274d-47e7-a655-6651ad867082
# ╟─b3df876e-bda6-4cc4-b305-f816168a4ae7
# ╠═b21feaa1-89c2-4ef0-8744-8c8d6c95147e
# ╟─a4972591-26bd-4555-b96b-765f20051763
# ╠═7e9229d8-08ff-4442-b028-684e52f4f668
# ╟─b01fa482-2d15-470c-9e1d-35bf12fb7609
# ╠═acebf498-8cdf-475e-9609-f60e68d4a355
# ╠═90522d32-25d7-45b7-aca5-a0c04da39432
# ╠═7c45e2c6-7d6a-42c2-afaf-4fd4d110b4bb
# ╠═943902fb-6853-46d7-8a0c-c4ef3bcd8229
# ╠═b3475222-d143-481f-85a5-82ee6272a6fa
# ╠═29603238-6c17-4668-8594-d240a76a23b9
# ╠═deaed5af-5700-418f-a6f1-05e0c0637d75

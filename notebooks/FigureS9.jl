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
# Figure S9. Domain SST and RH
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

# ╔═╡ 7c45e2c6-7d6a-42c2-afaf-4fd4d110b4bb
begin
	gds = NCDataset(datadir("wrf","grid.nc"))
	wln = gds["longitude"][:,:]
	wlt = gds["latitude"][:,:]
	close(gds)
end

# ╔═╡ 943902fb-6853-46d7-8a0c-c4ef3bcd8229
begin
	ds  = NCDataset(datadir("FigureS9data.nc"))
	sst = ds["SST"][:,:]
	rh  = ds["RH"][:,:]
	close(ds)
end

# ╔═╡ deaed5af-5700-418f-a6f1-05e0c0637d75
begin
	pplt.close(); fig,axs = pplt.subplots(ncols=2,axwidth=1.5)

	c1 = axs[1].pcolormesh(wln,wlt,sst,extend="both",levels=297:0.5:302)
	c2 = axs[2].pcolormesh(wln,wlt,rh,extend="both",cmap="drywet",levels=75:95)

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

	fig.savefig(projectdir("figures","figS9-sstrh.png"),transparent=false,dpi=400)
	load(projectdir("figures","figS9-sstrh.png"))
end

# ╔═╡ Cell order:
# ╟─fa2f8740-f813-11ec-00e1-112e2dfacda7
# ╟─1d7446ba-f464-44df-89e2-ae2a5726e849
# ╟─ab294fae-101f-4587-a2f4-7d72254dd421
# ╟─ccdf14a6-b200-42db-8455-0ea4eeb5ae2d
# ╟─eea4a60c-274d-47e7-a655-6651ad867082
# ╟─7c45e2c6-7d6a-42c2-afaf-4fd4d110b4bb
# ╟─943902fb-6853-46d7-8a0c-c4ef3bcd8229
# ╟─deaed5af-5700-418f-a6f1-05e0c0637d75

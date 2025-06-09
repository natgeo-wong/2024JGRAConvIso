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
	using ETOPO
	using GeoRegions
	using NCDatasets
	using PlutoUI
	using Printf

	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("ultraplot")

	include(srcdir("common.jl"))

	md"Loading modules for the 2024JGRAConvIso project..."
end

# ╔═╡ fa2f8740-f813-11ec-00e1-112e2dfacda7
md"
# Figure S5. Valid GeoRegions obtained from Moisture Budget Balancing

In this notebook, we define additional GeoRegions of interest for plotting and for analysis based on WRF modelling output and as necessary for figures.
"

# ╔═╡ ccdf14a6-b200-42db-8455-0ea4eeb5ae2d
TableOfContents()

# ╔═╡ a6ab92dd-1fae-4a04-b8f3-6782ea67c60b
md"
### A. Loading Basic Datasets and Station Information
"

# ╔═╡ 189e1048-c92d-457e-a30e-f4e523b80afc
begin
	infoall = stninfoall()
	infody  = stninfody(); ndystn = size(infody,1)
	infomo  = stninfomo()
	md"Loading station location information ..."
end

# ╔═╡ 32c650df-ccd2-4adf-a3b7-56611fff1b46
begin
	coast = readdlm(datadir("coast.cst"),comments=true,comment_char='#')
	x = coast[:,1]
	y = coast[:,2]
	md"Loading coastlines data ..."
end

# ╔═╡ 42b3b68f-a787-428a-8133-31b5de3723bb
etpd = ETOPODataset(path=datadir())

# ╔═╡ a4880c4b-62f5-482a-9639-70a16698e762
geo = GeoRegion("OTREC_wrf_d02",path=srcdir())

# ╔═╡ 1e3ad9d7-9d95-4fb5-a9b5-25c2c1a4b1c8
lsd = getLandSea(etpd,geo)

# ╔═╡ f6e3ea47-7316-48ec-9dbe-aceaa2059a5d
md"
### B. GeoRegions
"

# ╔═╡ 912a101b-7eb9-4322-a713-031aeffff20d
begin
	geovec = Array{GeoRegion}(undef,12,5)
	for istn = 1 : 12
		geovec[istn,1] = GeoRegion("OTREC_wrf_stn$(@sprintf("%02d",istn))",path=srcdir())
		for ibox = 1 : 4
			geovec[istn,ibox+1] = GeoRegion("OTREC_wrf_stn$(@sprintf("%02d",istn))_box$(@sprintf("%d",ibox))",path=srcdir())
		end
	end
end

# ╔═╡ 3aa94609-b566-457b-a7d2-71a9f2789825
begin
	wvc = readdlm(datadir("wrf","wrfvscalc-smooth30.txt"))
	qvl = readdlm(datadir("wrf","qbudgetdiff-smooth30.txt"))
end

# ╔═╡ cb306a1f-affd-432c-8c46-977758531654
begin
	pplt.close(); fig,axs = pplt.subplots(ncols=2,axwidth=2,sharey=0)

	c = 
	axs[1].pcolormesh(lsd.lon,lsd.lat,lsd.z'./1000,cmap="bukavu",levels=-2:0.2:2,extend="both")
	axs[2].pcolormesh(lsd.lon,lsd.lat,lsd.z'./1000,cmap="bukavu",levels=-2:0.2:2,extend="both")
	axs[1].plot([0,1],[0,1],c="r",lw=1,label="Valid Region",legend="t",legend_kw=Dict("frame"=>false,"ncols"=>1))
	axs[2].plot([0,1],[0,1],c="b",lw=1,linestyle="--",label="Invalid Region",legend="t",legend_kw=Dict("frame"=>false,"ncols"=>1))
	
	for ax in axs
		for istn = 1 : 12
			clon,clat = geovec[istn,1].geometry.centroid
			slon,slat = coordinates(geovec[istn,1])
			ax.scatter(clon,clat,c="r")
			(wvc[istn,1] < 0.15) .& (qvl[istn,1] < 0.05) ? ax.plot(slon,slat,c="r",lw=1) : ax.plot(slon,slat,c="b",lw=1,linestyle="--")
			for ibox = 1 : 4
				slon,slat = coordinates(geovec[istn,ibox+1])
				(wvc[istn,ibox+1] < 0.15) .& (qvl[istn,ibox+1] < 0.05) ? ax.plot(slon,slat,c="r",lw=1) : ax.plot(slon,slat,c="b",lw=1,linestyle="--")
			end
		end
		ax.format(xlabel=L"Longitude / $\degree$")
	end
	
	axs[1].format(xlim=(-80,-75),ylim=(2.5,7.5),ylabel=L"Latitude / $\degree$")
	axs[2].format(xlim=(-87,-82),ylim=(7,12))

	fig.colorbar(c,label="Topographic Height / km")
	fig.savefig(projectdir("figures","figS5-validboxes.png"),transparent=false,dpi=150)
	load(projectdir("figures","figS5-validboxes.png"))
end

# ╔═╡ Cell order:
# ╟─fa2f8740-f813-11ec-00e1-112e2dfacda7
# ╟─1d7446ba-f464-44df-89e2-ae2a5726e849
# ╟─ab294fae-101f-4587-a2f4-7d72254dd421
# ╟─ccdf14a6-b200-42db-8455-0ea4eeb5ae2d
# ╟─a6ab92dd-1fae-4a04-b8f3-6782ea67c60b
# ╟─189e1048-c92d-457e-a30e-f4e523b80afc
# ╟─32c650df-ccd2-4adf-a3b7-56611fff1b46
# ╟─42b3b68f-a787-428a-8133-31b5de3723bb
# ╟─a4880c4b-62f5-482a-9639-70a16698e762
# ╟─1e3ad9d7-9d95-4fb5-a9b5-25c2c1a4b1c8
# ╟─f6e3ea47-7316-48ec-9dbe-aceaa2059a5d
# ╟─912a101b-7eb9-4322-a713-031aeffff20d
# ╟─3aa94609-b566-457b-a7d2-71a9f2789825
# ╟─cb306a1f-affd-432c-8c46-977758531654

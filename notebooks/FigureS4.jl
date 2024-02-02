### A Pluto.jl notebook ###
# v0.19.37

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

	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")

	include(srcdir("common.jl"))

	md"Loading modules for the 2024GLConvIso project..."
end

# ╔═╡ fa2f8740-f813-11ec-00e1-112e2dfacda7
md"
# Figure 2. OTREC Region of Interest
"

# ╔═╡ ccdf14a6-b200-42db-8455-0ea4eeb5ae2d
TableOfContents()

# ╔═╡ a6ab92dd-1fae-4a04-b8f3-6782ea67c60b
md"
### A. Loading Datasets and LandSea Masks ...
"

# ╔═╡ 189e1048-c92d-457e-a30e-f4e523b80afc
begin
	infoall = stninfoall()
	infody  = stninfody()
	infomo  = stninfomo()
	infocr  = stninfocostarica()
	md"Loading station location information ..."
end

# ╔═╡ 32c650df-ccd2-4adf-a3b7-56611fff1b46
begin
	coast = readdlm(datadir("coast.cst"),comments=true)
	x = coast[:,1]
	y = coast[:,2]
	md"Preloading coastline data"
end

# ╔═╡ eb7a010b-5e93-48c5-8e0a-698fbd70cda4
begin
	flights = loadflightpaths()
	md"Loading flight path data ..."
end

# ╔═╡ 60d38ccd-6538-426f-985e-9d782834dce9
begin
	xinset = deepcopy(x); xinset[(x.>285).|(x.<270).|(y.>15).|(y.<0)].=NaN
	yinset = deepcopy(y); yinset[(x.>285).|(x.<270).|(y.>15).|(y.<0)].=NaN
	md"Filtering out coastlines outside domain ..."
end

# ╔═╡ c16d89d9-d7ba-4b79-aa7c-8570467333e0
geo = GeoRegion("OTREC")

# ╔═╡ 5d6c3cd6-e406-461d-a226-20022060398d
lsd = getLandSea(geo,path=datadir(),returnlsd=true)

# ╔═╡ bc8992d5-fea2-405d-bf6d-9e83c0abf768
tgeo = RectRegion("TMP","OTREC","TMP",[7,3,285,281],savegeo=false)

# ╔═╡ 982993ab-3437-4e05-bb8c-a78834d4c39d
ggrd = RegionGrid(tgeo,lsd.lon,lsd.lat)

# ╔═╡ ff2ded99-0764-47d5-a2cc-36d5b56acbc2
begin
	z = extractGrid(lsd.z,ggrd)
	md"Extracting data for temporary region"
end

# ╔═╡ 1a643e58-39c1-4c6b-b340-978056871b6b
md"
### C. Plotting Regions of Interest
"

# ╔═╡ 6dfab570-8d7b-47e2-9971-69fde82cbe70
begin
	blon_02,blat_02 = coordGeoRegion(GeoRegion("OTREC_STN02"))
	blon_03,blat_03 = coordGeoRegion(GeoRegion("OTREC_STN03"))
	blon_04,blat_04 = coordGeoRegion(GeoRegion("OTREC_STN04"))
	blon_02_01,blat_02_01 = coordGeoRegion(GeoRegion("OTREC_STN02_01"))
	blon_03_01,blat_03_01 = coordGeoRegion(GeoRegion("OTREC_STN03_01"))
	blon_04_01,blat_04_01 = coordGeoRegion(GeoRegion("OTREC_STN04_01"))
	blon_02_04,blat_02_04 = coordGeoRegion(GeoRegion("OTREC_STN02_04"))
	blon_03_04,blat_03_04 = coordGeoRegion(GeoRegion("OTREC_STN03_04"))
	blon_04_04,blat_04_04 = coordGeoRegion(GeoRegion("OTREC_STN04_04"))
	md"Loading GeoRegion bound coordinates ..."
end

# ╔═╡ 940fd96a-98cf-441e-bbba-20ee1250551d
clrs = pplt.get_colors("default")

# ╔═╡ d7755534-3565-4011-b6e3-e131991008db
begin
	pplt.close(); fig,axs = pplt.subplots(aspect=1,axwidth=3)
	textdict = Dict("fc"=>"grey3","ec"=>"none","alpha"=>0.6)

	c = axs[1].pcolormesh(
		ggrd.lon,ggrd.lat,z'/1000,
		levels=-4.5:0.5:4.5,cmap="delta",extend="both"
	)
	axs[1].contour(ggrd.lon,ggrd.lat,z'/1000,levels=[0],color="k",lw=1)
	axs[1].scatter(infody[:,2],infody[:,3],c="r",zorder=5)
	axs[1].plot(blon_02,blat_02,c=clrs[1])
	axs[1].plot(blon_03,blat_03,c=clrs[2])
	axs[1].plot(blon_04,blat_04,c=clrs[3])
	axs[1].plot(blon_02_01,blat_02_01,c=clrs[1],linestyle=":")
	axs[1].plot(blon_03_01,blat_03_01,c=clrs[2],linestyle=":")
	axs[1].plot(blon_04_01,blat_04_01,c=clrs[3],linestyle=":")
	axs[1].plot(blon_02_04,blat_02_04,c=clrs[1],linestyle=":")
	axs[1].plot(blon_03_04,blat_03_04,c=clrs[2],linestyle=":")
	axs[1].plot(blon_04_04,blat_04_04,c=clrs[3],linestyle=":")
	axs[1].text(281.4,5.9,"Bahía Solano",c="k",size=8,bbox=textdict)
	axs[1].text(283.6,5.4,"Quibdó",c="k",size=8,bbox=textdict)
	axs[1].text(283.3,3.5,"Buenaventura",c="k",size=8,bbox=textdict)
	axs[1].format(xlim=(281,285),ylim=(3,7))

	fig.colorbar(c,label="Topographic Height/ km")
	fig.savefig(projectdir("figures","figS4-georegionbreakdown.png"),transparent=false,dpi=400)
	load(projectdir("figures","figS4-georegionbreakdown.png"))
end

# ╔═╡ Cell order:
# ╟─fa2f8740-f813-11ec-00e1-112e2dfacda7
# ╟─1d7446ba-f464-44df-89e2-ae2a5726e849
# ╟─ab294fae-101f-4587-a2f4-7d72254dd421
# ╟─ccdf14a6-b200-42db-8455-0ea4eeb5ae2d
# ╟─a6ab92dd-1fae-4a04-b8f3-6782ea67c60b
# ╟─189e1048-c92d-457e-a30e-f4e523b80afc
# ╟─32c650df-ccd2-4adf-a3b7-56611fff1b46
# ╟─eb7a010b-5e93-48c5-8e0a-698fbd70cda4
# ╟─60d38ccd-6538-426f-985e-9d782834dce9
# ╟─c16d89d9-d7ba-4b79-aa7c-8570467333e0
# ╟─5d6c3cd6-e406-461d-a226-20022060398d
# ╟─bc8992d5-fea2-405d-bf6d-9e83c0abf768
# ╟─982993ab-3437-4e05-bb8c-a78834d4c39d
# ╟─ff2ded99-0764-47d5-a2cc-36d5b56acbc2
# ╟─1a643e58-39c1-4c6b-b340-978056871b6b
# ╟─6dfab570-8d7b-47e2-9971-69fde82cbe70
# ╠═940fd96a-98cf-441e-bbba-20ee1250551d
# ╟─d7755534-3565-4011-b6e3-e131991008db

### A Pluto.jl notebook ###
# v0.20.4

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
	using ETOPO
	using GeoRegions
	using NCDatasets
	using PlutoUI
	using Printf

	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("ultraplot")

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

# ╔═╡ eb7a010b-5e93-48c5-8e0a-698fbd70cda4
begin
	flights = loadflightpaths()
	md"Loading flight path data ..."
end

# ╔═╡ 3ec67c8a-bbbf-4ad5-a55e-e681c88d3c9d
etpd = ETOPODataset(path=datadir())

# ╔═╡ d5adfc26-72a1-49c0-b8bf-f626bd4e69bc
geo_d01 = GeoRegion("OTREC_wrf_d01",path=srcdir())

# ╔═╡ c16d89d9-d7ba-4b79-aa7c-8570467333e0
geo_d02 = GeoRegion("OTREC_wrf_d02",path=srcdir())

# ╔═╡ 32c650df-ccd2-4adf-a3b7-56611fff1b46
begin
	coast = readdlm(datadir("coast.cst"),comments=true)
	x = coast[:,1]; y = coast[:,2]; npnts = length(x)
	for ipnt = 1 : npnts
		pnt = Point2(x[ipnt],y[ipnt])
		if !in(pnt,geo_d02)
			x[ipnt] = NaN; y[ipnt] = NaN
		end
	end
	md"Preloading coastline data"
end

# ╔═╡ 194bb86d-a3bb-4585-ba66-067dd4488189
geo_big = GeoRegion("Fig2_big",path=srcdir())

# ╔═╡ ee60cad1-0863-4974-9690-087d0c6dc06e
geo_sml = GeoRegion("Fig2_small",path=srcdir())

# ╔═╡ 1ae33598-e283-49f7-bc8d-274b33ab12a5
lsd_d02 = getLandSea(etpd,geo_d02,returnlsd=true)

# ╔═╡ 5d6c3cd6-e406-461d-a226-20022060398d
lsd_sml = getLandSea(etpd,geo_sml,returnlsd=true)

# ╔═╡ 76e35851-71a0-4b89-98bb-9a82bf34bd34
lsd_big = getLandSea(etpd,geo_big,returnlsd=true)

# ╔═╡ 7b4eaccc-7f7e-4c82-9ec5-f82834c75098
lon_d01,lat_d01 = coordinates(geo_d01);

# ╔═╡ dc660494-dc7b-41bc-8ca8-ae68c2bfd94b
lon_d02,lat_d02 = coordinates(geo_d02);

# ╔═╡ 1a643e58-39c1-4c6b-b340-978056871b6b
md"
### C. Plotting Regions of Interest
"

# ╔═╡ d7755534-3565-4011-b6e3-e131991008db
begin
	pplt.close(); fig,axs = pplt.subplots(aspect=27/17,axwidth=5)

	lvls = -6 : 6
	textdict = Dict("fc"=>"grey3","ec"=>"none","alpha"=>0.6)
	
	axs[1].pcolormesh(
		lsd_big.lon[1:10:end],lsd_big.lat[1:10:end],
		lsd_big.z[1:10:end,1:10:end]'/1000,
		alpha=0.3,levels=lvls,cmap="delta",extend="both"
	)
	
	axs[1].pcolormesh(
		lsd_d02.lon.+360,lsd_d02.lat,lsd_d02.z'/1000,
		levels=lvls,cmap="delta",extend="both"
	)

	for ipnt = 1 : 13
		geosample = GeoRegion("OTREC_wrf_ITCZ$(@sprintf("%02d",ipnt))",path=srcdir())
		ilon,ilat = coordinates(geosample)
		axs[1].plot(ilon.+360,ilat,c="yellow")
	end
	for ipnt = 1 : 13
		geosample = GeoRegion("OTREC_wrf_PAC2ATL$(@sprintf("%02d",ipnt))",path=srcdir())
		ilon,ilat = coordinates(geosample)
		axs[1].plot(ilon.+360,ilat,c="brown")
	end
	for ipnt = 1 : 11
		geosample = GeoRegion("OTREC_wrf_CrossITCZ$(@sprintf("%02d",ipnt))",path=srcdir())
		ilon,ilat = coordinates(geosample)
		axs[1].plot(ilon.+360,ilat,c="b")
	end
	
	axs[1].plot(x,y,lw=0.5,c="k")
	axs[1].plot([275,275,278,278,275],[8,10.5,10.5,8,8],lw=1,c="k",linestyle="--")
	axs[1].scatter(infody[:,2],infody[:,3],zorder=4,c="pink")
	axs[1].scatter(infocr[:,2],infocr[:,3],zorder=4,c="red")
	axs[1].plot([275,288],[10.5,15],lw=1,c="k",linestyle=":")
	axs[1].plot([278,295],[8,9.03],lw=1,c="k",linestyle=":")
	axs[1].plot(lon_d02.+360,lat_d02,lw=5,c="k")

	axs[1].text(272.3,10,"Liberia",c="k",size=8,bbox=textdict)
	axs[1].text(277.5,13,"San Andres",c="k",size=8,bbox=textdict)
	axs[1].text(278.6,6,"Bahía Solano",c="k",size=8,bbox=textdict)
	axs[1].text(282,4.75,"Quibdó",c="k",size=8,bbox=textdict)
	axs[1].text(278.8,3.5,"Buenaventura",c="k",size=8,bbox=textdict)
	
	axs[1].format(
		xlim=(269,296),xlocator=255:5:300,xminorlocator=240:2.5:315,
		ylim=(-1,16),ylocator=-15:5:30,yminorlocator=-20:2.5:35,
		xlabel=L"Longitude / $\degree$",ylabel=L"Latitude / $\degree$",
		grid=true,gridcolor="w",#suptitle="Available Colombia Stations",
	)

	ix = fig.add_axes([0.635,0.632,0.21,0.30])
	c = ix.pcolormesh(
		lsd_sml.lon,lsd_sml.lat,lsd_sml.z'/1000,
		levels=lvls,cmap="delta",extend="both"
	)
	ix.plot(x,y,lw=0.5,c="k")
	ix.scatter(infocr[:,2],infocr[:,3],s=10,zorder=4,c="r")

	ix.text(275.2,10.15,"EEFMB",c="k",size=7,bbox=textdict)
	ix.text(275.4,9.65,"CGFI",c="k",size=7,bbox=textdict)
	ix.text(276.1,9.64,"AMDQ",c="k",size=7,bbox=textdict)
	ix.text(276.1,8.3,"OSA",c="k",size=7,bbox=textdict)
	ix.text(276.2,10.22,"Bataan",c="k",size=7,bbox=textdict)
	ix.text(277,10.1,"Limon",c="k",size=7,bbox=textdict)
	ix.text(277,9.45,"Cahuita",c="k",size=7,bbox=textdict)
	ix.text(275.1,8.2,"(a)",c="k",size=9)
	
	ix.format(
		xlim=(275,278),ylim=(8,10.5),xtickloc="none",ytickloc="none",
		xlocator=255:15:300,xminorlocator=240:5:315,xticklabels=[],
		ylocator=-15:15:30,yminorlocator=-30:5:45,yticklabels=[]
	)

	ix = fig.add_axes([0.635,0.175,0.21,0.36])
	ix.pcolormesh(
		lsd_big.lon[1:10:end],lsd_big.lat[1:10:end],
		lsd_big.z[1:10:end,1:10:end]'/1000,alpha=0.3,
		levels=lvls,cmap="delta",extend="both"
	)
	ix.pcolormesh(
		lsd_d02.lon.+360,lsd_d02.lat,lsd_d02.z'/1000,
		levels=lvls,cmap="delta",extend="both"
	)

	ix.plot(x,y,lw=0.5,c="k")
	ix.plot(lon_d02.+360,lat_d02,lw=1.5,c="k")
	ix.plot(lon_d01.+360,lat_d01,lw=1.5,c="k",linestyle="--")
	ix.text(281,-5,"d02",c="k",size=8)
	ix.text(286,-14,"d01",c="k",size=8)
	ix.text(252,28,"(b)",c="k",size=9)
	ix.format(
		xlim=(250,305),ylim=(-20,35),xtickloc="none",ytickloc="none",
		xlocator=255:15:300,xminorlocator=240:5:315,xticklabels=[],
		ylocator=-15:15:30,yminorlocator=-30:5:45,yticklabels=[]
	)

	fig.colorbar(c,label="Topographic Height/ km",length=0.8)
	fig.savefig(projectdir("figures","fig2-domain.png"),transparent=false,dpi=400)
	load(projectdir("figures","fig2-domain.png"))
end

# ╔═╡ Cell order:
# ╟─fa2f8740-f813-11ec-00e1-112e2dfacda7
# ╟─1d7446ba-f464-44df-89e2-ae2a5726e849
# ╟─ab294fae-101f-4587-a2f4-7d72254dd421
# ╟─ccdf14a6-b200-42db-8455-0ea4eeb5ae2d
# ╟─a6ab92dd-1fae-4a04-b8f3-6782ea67c60b
# ╟─189e1048-c92d-457e-a30e-f4e523b80afc
# ╠═32c650df-ccd2-4adf-a3b7-56611fff1b46
# ╟─eb7a010b-5e93-48c5-8e0a-698fbd70cda4
# ╠═3ec67c8a-bbbf-4ad5-a55e-e681c88d3c9d
# ╟─d5adfc26-72a1-49c0-b8bf-f626bd4e69bc
# ╟─c16d89d9-d7ba-4b79-aa7c-8570467333e0
# ╟─194bb86d-a3bb-4585-ba66-067dd4488189
# ╟─ee60cad1-0863-4974-9690-087d0c6dc06e
# ╟─1ae33598-e283-49f7-bc8d-274b33ab12a5
# ╟─5d6c3cd6-e406-461d-a226-20022060398d
# ╠═76e35851-71a0-4b89-98bb-9a82bf34bd34
# ╠═7b4eaccc-7f7e-4c82-9ec5-f82834c75098
# ╠═dc660494-dc7b-41bc-8ca8-ae68c2bfd94b
# ╟─1a643e58-39c1-4c6b-b340-978056871b6b
# ╟─d7755534-3565-4011-b6e3-e131991008db

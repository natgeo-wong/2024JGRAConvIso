### A Pluto.jl notebook ###
# v0.19.37

using Markdown
using InteractiveUtils

# ╔═╡ 9db63111-71ec-43b9-9174-d4921191eafd
begin
	using Pkg; Pkg.activate()
	using DrWatson
	md"Using DrWatson to ensure reproducibility between different machines ..."
end

# ╔═╡ da16ef31-2091-4c36-b3d0-05f4197771f6
begin
	@quickactivate "2024GLConvIso"
	using DelimitedFiles
	using ERA5Reanalysis
	using NCDatasets
	using Printf
	using Statistics
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")

	include(srcdir("common.jl"))
	
	md"Loading modules for the 2024GLConvIso project..."
end

# ╔═╡ e0a4eb5e-467e-11ec-21da-af325ae0fd15
md"
# 02c. GPM Climatology for various OTREC Stations

Tex
"

# ╔═╡ 19b04d89-520c-4fdc-8152-3e410ae08b8f
begin
	coast = readdlm(datadir("coast.cst"),comments=true)
	clon  = coast[:,1]
	clat  = coast[:,2]
	md"Preloading coastline data"
end

# ╔═╡ a567368e-ec0e-4944-9c1f-be1cfe9c5d1d
md"
### A. Loading $p_\omega$ Data
"

# ╔═╡ 4c887a8c-c057-4695-8126-8881a92220a8
begin
	eds = NCDataset(datadir("p_wwgt-compiledera5.nc"))
	lon = eds["longitude"][:]
	lat = eds["latitude"][:]
	pωe = eds["p_wwgt"][:,:]
	close(eds)
	md"Loading ERA5 $p_\omega$ data ..."
end

# ╔═╡ 0cbc9af9-5bbc-484d-b9c0-0474831c5ced
begin
	wds = NCDataset(datadir("wrf","grid.nc"))
	wln = wds["longitude"][:,:]
	wlt = wds["latitude"][:,:]
	close(wds)
	wds = NCDataset(datadir("p_wwgt-compiledwrf.nc"))
	pωw  = wds["p_wwgt"][:,:]
	close(wds)
	md"Loading WRF $p_\omega$ data ..."
end

# ╔═╡ a443e74d-263b-4c83-8cbd-f482cda0ff79
md"
### B. Extracting ERA5 Tropical Data into OTREC Region
"

# ╔═╡ 4de46224-f340-4437-b4c8-7b53440f74d0
e5ds = ERA5Dummy(path=datadir())

# ╔═╡ 9d55b585-71d0-42ff-b045-5f7ddde29ba6
geo = GeoRegion("OTREC")

# ╔═╡ b06ef4a2-c61d-4fb4-96f3-aaf9ceb1d8a4
egeo = ERA5Region(geo)

# ╔═╡ a92c0cdd-fcea-4bcb-ac1d-e467eb26da6b
lsd = getLandSea(e5ds,egeo)

# ╔═╡ 336a4b18-bc1a-4c76-88e0-578f9971c5da
ggrd = RegionGrid(geo,lon,lat)

# ╔═╡ ef6d666b-8a64-461e-9fb5-440748135683
begin
	pω_OTREC = extractGrid(pωe,ggrd)
	pω_OTREC[lsd.z.>500] .= NaN
	md"Extracting ERA5 $p_\omega$ for OTREC Region"
end

# ╔═╡ 7e69c41b-6e5e-4d5b-8cda-2b5d65ddf6c5
md"
### C. Binning WRF into OTREC ERA5 Grid
"

# ╔═╡ 93fa8a76-03d6-47b7-8e01-6c72f3af09de
begin
	ipnt_lon = zeros(Int,543,543)
	ipnt_lat = zeros(Int,543,543)
	ind      = zeros(Bool,543,543)
	for ilat = 1 : 543, ilon = 1 : 543
		ipnt_lon[ilon,ilat] = argmin(abs.(wln[ilon,ilat].-ggrd.lon.+360))
		ipnt_lat[ilon,ilat] = argmin(abs.(wlt[ilon,ilat].-ggrd.lat))
	end
	md"Finding closest ERA5 points to each of the WRF points ..."
end

# ╔═╡ f2d6e36b-e12f-4c24-81a8-5e03395add53
begin
	wwgtp = zeros(length(ggrd.lon),length(ggrd.lat))
	md"Preallocation of arrays for weighted pressure ..."
end

# ╔═╡ 0972ac7d-235b-4121-9e40-210aac89ecad
begin
	for ilat = 1 : length(ggrd.lat), ilon = 1 : length(ggrd.lon)
		ind = (ipnt_lon.==ilon).&(ipnt_lat.==ilat)
		iipωw = @view pωw[ind]
		wwgtp[ilon,ilat] = mean(iipωw[.!isnan.(iipωw)])
	end
	wwgtp[lsd.z.>500] .= NaN
	md"Regridding/binning $p_\omega$ into ERA5 Grid"
end

# ╔═╡ 4d5f7822-4ccf-406f-84d7-1edd8cd3d38b
md"
### D. Plotting
"

# ╔═╡ b87f043f-3994-4587-bfde-92f2808d27e7
begin
	pplt.close(); fig,axs = pplt.subplots(ncols=3,axwidth=1.5)
	
	c1 = axs[1].pcolormesh(
		ggrd.lon,ggrd.lat,pω_OTREC'/100,
		levels=450:25:750,cmap="drywet_r",extend="both"
	)
	axs[2].pcolormesh(
		ggrd.lon,ggrd.lat,wwgtp'/100,
		levels=450:25:750,cmap="drywet_r",extend="both"
	)
	c2 = axs[3].pcolormesh(
		ggrd.lon,ggrd.lat,wwgtp'/100 .- pω_OTREC'/100,
		levels=-150:25:150,cmap="rdbu",extend="both"
	)

	axs[1].format(ltitle="(a) ERA5")
	axs[2].format(ltitle="(b) WRF")
	axs[3].format(ltitle="(c) WRF - ERA5")

	for ax in axs
		ax.plot(clon,clat,c="k",lw=1)
		ax.format(
			xlim=(270,285),xlocator=270:5:285,xlabel=L"Longitude / $\degree$",
			ylim=(0,15),ylocator=0:5:15,ylabel=L"Latitude / $\degree$"
		)
	end

	axs[2].colorbar(c1,label=L"$p_\omega$ / hPa")
	axs[3].colorbar(c2,label=L"$\Delta p_\omega$ / hPa")
	
	fig.savefig(projectdir("figures","figS3-era5vswrf.png"),transparent=false,dpi=400)
	load(projectdir("figures","figS3-era5vswrf.png"))
end

# ╔═╡ Cell order:
# ╟─e0a4eb5e-467e-11ec-21da-af325ae0fd15
# ╟─9db63111-71ec-43b9-9174-d4921191eafd
# ╟─da16ef31-2091-4c36-b3d0-05f4197771f6
# ╟─19b04d89-520c-4fdc-8152-3e410ae08b8f
# ╟─a567368e-ec0e-4944-9c1f-be1cfe9c5d1d
# ╟─4c887a8c-c057-4695-8126-8881a92220a8
# ╟─0cbc9af9-5bbc-484d-b9c0-0474831c5ced
# ╟─a443e74d-263b-4c83-8cbd-f482cda0ff79
# ╟─4de46224-f340-4437-b4c8-7b53440f74d0
# ╟─9d55b585-71d0-42ff-b045-5f7ddde29ba6
# ╟─b06ef4a2-c61d-4fb4-96f3-aaf9ceb1d8a4
# ╟─a92c0cdd-fcea-4bcb-ac1d-e467eb26da6b
# ╟─336a4b18-bc1a-4c76-88e0-578f9971c5da
# ╟─ef6d666b-8a64-461e-9fb5-440748135683
# ╟─7e69c41b-6e5e-4d5b-8cda-2b5d65ddf6c5
# ╟─93fa8a76-03d6-47b7-8e01-6c72f3af09de
# ╟─f2d6e36b-e12f-4c24-81a8-5e03395add53
# ╟─0972ac7d-235b-4121-9e40-210aac89ecad
# ╟─4d5f7822-4ccf-406f-84d7-1edd8cd3d38b
# ╟─b87f043f-3994-4587-bfde-92f2808d27e7

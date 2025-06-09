### A Pluto.jl notebook ###
# v0.20.10

using Markdown
using InteractiveUtils

# ╔═╡ b3bc56eb-43a7-4736-bd66-704529911d60
begin
	using Pkg; Pkg.activate()
	using DrWatson
	md"Using DrWatson to ensure reproducibility between different machines ..."
end

# ╔═╡ 9071763e-f6ad-4468-ae48-369307a85263
begin
	@quickactivate "2024GLConvIso"
	using Dates
	using DelimitedFiles
	using ERA5Reanalysis
	using NASAPrecipitation
	using NCDatasets
	using Printf
	using Statistics
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("ultraplot")

	include(srcdir("common.jl"))
	
	md"Loading modules for the 2024GLConvIso project..."
end

# ╔═╡ a8431d46-46fe-11ec-2b8d-e39caffdabec
md"
# Figure 3a. Rainfall, IMERGv7 vs ERA5 vs WRF

In this notebook, we explicitly compare the rainfall measured by the 3 stations we have, to GPM Rainfall data at the nearest gridpoint corresponding to these stations.
"

# ╔═╡ 4579f773-e2b4-4b1f-a8c1-485341d451fc
begin
	coast = readdlm(datadir("coast.cst"),comments=true,comment_char='#')
	x = coast[:,1]
	y = coast[:,2]
	md"Loading coastlines data ..."
end

# ╔═╡ c0799dcd-7c88-4f75-ac19-76b615ff4454
e5ds = ERA5Hourly(start=Date(2019,8,1),stop=Date(2020,12,31),path=datadir())

# ╔═╡ 9a455f45-592c-43ab-a71c-d1de27bf8552
egeo = ERA5Region("OTREC_wrf_d02",path=srcdir())

# ╔═╡ 5260b50e-7436-4afa-96b4-3e9ac1d26229
elsd = getLandSea(e5ds,egeo)

# ╔═╡ 6c91263c-a1e2-41f6-83f1-e69e55ef0c38
begin
	wds = NCDataset(datadir("wrf","processed","p_wwgt-compiledwrf-20190801_20201231.nc"))
	wpw = wds["p_wwgt"][:,:] ./ 100
	wrn = wds["RAINNC"][:,:]
	close(wds)
end

# ╔═╡ 4b0ea8d0-4588-4abf-88cb-4dfd8e52b4a7
begin
	gds = NCDataset(datadir("wrf","grid.nc"))
	wln = gds["longitude"][:,:]; nx = size(wln,1)
	wlt = gds["latitude"][:,:];  ny = size(wlt,2)
	close(gds)
end

# ╔═╡ d68d2d48-7f33-44bb-8d40-996167585e45
begin
	eds = NCDataset(datadir("p_wwgt-compiledera5.nc"))
	eln = eds["longitude"][:]
	elt = eds["latitude"][:]
	epw = eds["p_wwgt"][:,:] ./ 100
	ggrd = RegionGrid(egeo,eln,elt)
	epw = extract(epw,ggrd)
	close(eds)
end

# ╔═╡ 59e0c4ad-c89a-4d70-9aa2-72a5c3da720e
begin
	wepw  = zeros(length(elsd.lon),length(elsd.lat))
	wern  = zeros(length(elsd.lon),length(elsd.lat))
	ipnt_lon = zeros(Int,nx,ny)
	ipnt_lat = zeros(Int,nx,ny)
	for ilat = 1 : ny, ilon = 1 : nx
		ipnt_lon[ilon,ilat] = argmin(abs.(wln[ilon,ilat].-elsd.lon))
		ipnt_lat[ilon,ilat] = argmin(abs.(wlt[ilon,ilat].-elsd.lat))
	end
	for ilat = 1 : length(elsd.lat), ilon = 1 : length(elsd.lon)
		ind = (ipnt_lon.==ilon).&(ipnt_lat.==ilat)
		idata = @view wpw[ind]
		wepw[ilon,ilat] = mean(idata[.!isnan.(idata)])
		idata = @view wrn[ind]
		wern[ilon,ilat] = mean(idata[.!isnan.(idata)])
	end
	wepw[wern.<5] .= NaN
end

# ╔═╡ b74efe40-6288-4b27-b757-5d0771f2552e
begin
	pplt.close(); fig,axs = pplt.subplots(ncols=3,axwidth=1.5,aspect=1)
	
	c1 = 
	axs[1].pcolormesh(elsd.lon,elsd.lat,epw',levels=450:25:850,cmap="drywet_r",extend="both")
	axs[2].pcolormesh(elsd.lon,elsd.lat,wepw',levels=450:25:850,cmap="drywet_r",extend="both")
	c2 = 
	axs[3].pcolormesh(elsd.lon,elsd.lat,(wepw.-epw)',levels=-0:25:250,cmap="Reds",extend="both")

	axs[1].format(ltitle="(a) ERA5")
	axs[2].format(ltitle="(b) Regridded WRF")
	axs[3].format(ltitle="(c) WRF - ERA5")

	for ax in axs
		ax.plot(x,y,lw=0.5,c="k")
		ax.format(
			xlim=(-90,-75),xlabel=L"Longitude / $\degree$",
			ylim=(0,15),ylabel=L"Latitude / $\degree$",ylocator=0:5:15,
			# suptitle = "(a) Mean Rainfall Rate (2019 Aug - 2020 Dec)"
		)
	end

	axs[2].colorbar(c1,label=L"$p_\omega$ / hPa")
	axs[3].colorbar(c2,label=L"$\Delta p_\omega$ / hPa")
	
	fig.savefig(projectdir("figures","figS3-era5vswrf.png"),transparent=false,dpi=400)
	load(projectdir("figures","figS3-era5vswrf.png"))
end

# ╔═╡ Cell order:
# ╟─a8431d46-46fe-11ec-2b8d-e39caffdabec
# ╟─b3bc56eb-43a7-4736-bd66-704529911d60
# ╟─9071763e-f6ad-4468-ae48-369307a85263
# ╟─4579f773-e2b4-4b1f-a8c1-485341d451fc
# ╟─c0799dcd-7c88-4f75-ac19-76b615ff4454
# ╠═9a455f45-592c-43ab-a71c-d1de27bf8552
# ╟─5260b50e-7436-4afa-96b4-3e9ac1d26229
# ╠═6c91263c-a1e2-41f6-83f1-e69e55ef0c38
# ╠═4b0ea8d0-4588-4abf-88cb-4dfd8e52b4a7
# ╠═d68d2d48-7f33-44bb-8d40-996167585e45
# ╟─59e0c4ad-c89a-4d70-9aa2-72a5c3da720e
# ╟─b74efe40-6288-4b27-b757-5d0771f2552e

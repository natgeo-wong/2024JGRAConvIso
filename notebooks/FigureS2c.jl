### A Pluto.jl notebook ###
# v0.20.4

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
# Figure 2. Station Rainfall and Isotopes

In this notebook, we explicitly compare the rainfall measured by the 3 stations we have, to GPM Rainfall data at the nearest gridpoint corresponding to these stations.
"

# ╔═╡ 4579f773-e2b4-4b1f-a8c1-485341d451fc
begin
	coast = readdlm(datadir("coast.cst"),comments=true,comment_char='#')
	x = coast[:,1]
	y = coast[:,2]
	md"Loading coastlines data ..."
end

# ╔═╡ 7592bde2-b5dc-4b1e-8d51-4a9c33ba25d3
npd = IMERGMonthly(start=Date(2019),stop=Date(2020),path=datadir())

# ╔═╡ c0799dcd-7c88-4f75-ac19-76b615ff4454
e5ds = ERA5Monthly(start=Date(2019),stop=Date(2020),path=datadir())

# ╔═╡ dc809308-c782-4690-bd79-b4d0d2683cf1
evar = SingleVariable("p_wwgt",path=srcdir())

# ╔═╡ 538d96e5-ed7b-477a-8f3a-959a9b978ec5
geo = GeoRegion("OTREC_wrf_d02",path=srcdir())

# ╔═╡ 9a455f45-592c-43ab-a71c-d1de27bf8552
egeo = ERA5Region(geo)

# ╔═╡ a822f787-214a-4ef3-ba54-b899e890aeff
nlsd = getLandSea(npd,geo)

# ╔═╡ 5260b50e-7436-4afa-96b4-3e9ac1d26229
elsd = getLandSea(e5ds,egeo)

# ╔═╡ 9520065b-987c-4b08-8218-3aca80cd4d20
yr = 2020

# ╔═╡ 696d9f79-c66c-4ff7-9564-aa2fae5cf76b
mo = 12

# ╔═╡ cf67be7e-e020-4a47-8ddc-320f95194c7c
dtvec = Date(2019,8) : Month(1) : Date(2020,12)

# ╔═╡ d779d4b0-0c3b-4a56-bc44-0ada500c529b
begin
	ds   = NCDataset(datadir("wrf","grid.nc"))
    wlon = ds["longitude"][:,:]
    wlat = ds["latitude"][:,:]
    close(ds)
end

# ╔═╡ 19cd09a4-9f58-49cc-ba13-326ce1392413
begin
	wds = NCDataset(datadir("wrf","processed","p_wwgt-compiledwrf-monthly-20190801_20201231.nc"))
	wpwgt = wds["p_wwgt"][:,:,dtvec.==Date(yr,mo)] / 100
	wprcp = wds["RAINNC"][:,:,dtvec.==Date(yr,mo)]
	close(wds)
end

# ╔═╡ 42391f74-64d6-4972-9074-7e923e66c077
begin
	eds = read(e5ds,evar,egeo,Date(yr,mo))
	epwgt = eds[evar.ID][:,:,mo] / 100
	close(eds)
	eds = read(e5ds,SingleVariable("tp"),egeo,Date(yr,mo))
	eprcp = eds["tp"][:,:,mo] * 1000
	close(eds)
	epwgt[eprcp.<2.5] .= NaN
end

# ╔═╡ 1b46ee2b-3dfd-460f-b862-000ba7b51617
function regrid(wprcp,lsd)
	ds   = NCDataset(datadir("wrf","grid.nc"))
    wlon = ds["longitude"][:,:]; nx,ny = size(wlon)
    wlat = ds["latitude"][:,:]
    close(ds)
	ipnt_lon = zeros(Int,nx,ny)
	ipnt_lat = zeros(Int,nx,ny)
	for ilat = 1 : ny, ilon = 1 : nx
		ipnt_lon[ilon,ilat] = argmin(abs.(wlon[ilon,ilat].-lsd.lon))
		ipnt_lat[ilon,ilat] = argmin(abs.(wlat[ilon,ilat].-lsd.lat))
	end
	nlon = length(lsd.lon)
	nlat = length(lsd.lat)
	newprcp = zeros(nlon,nlat)
	for ilat = 1 : nlat, ilon = 1 : nlon
		ind = (ipnt_lon.==ilon).&(ipnt_lat.==ilat)
		iprcp = @view wprcp[ind]
		newprcp[ilon,ilat] = mean(iprcp[.!isnan.(iprcp)])
	end

	return newprcp
end

# ╔═╡ 330e270e-4203-4bfc-be20-98f3d9ae89fa
begin
	npwgt = regrid(wpwgt,elsd);
	nprcp = regrid(wprcp,elsd);
	npwgt[nprcp.<2.5] .= NaN
end

# ╔═╡ b74efe40-6288-4b27-b757-5d0771f2552e
begin
	asp = (geo.N-geo.S+2)/(geo.E-geo.W+2)
	pplt.close(); fig,axs = pplt.subplots(ncols=3,axwidth=1.5,aspect=asp)
	
	c1 = axs[1].pcolormesh(elsd.lon,elsd.lat,epwgt',levels=400:50:900,cmap="drywet_r",extend="both")
	axs[2].pcolormesh(elsd.lon,elsd.lat,npwgt',levels=400:50:900,cmap="drywet_r",extend="both")
	c2 = axs[3].pcolormesh(elsd.lon,elsd.lat,(npwgt.-epwgt)',levels=-200:50:200,cmap="RdBu_r",extend="both")

	axs[1].format(ultitle="(a) ERA5")
	axs[2].format(ultitle="(b) Regridded WRF")
	axs[3].format(ultitle="(c) WRF - ERA5")

	for ax in axs
		ax.plot(x,y,lw=0.5,c="k")
		ax.format(
			xlim=(geo.W-1,geo.E+1),xlabel=L"Longitude / $\degree$",
			ylim=(geo.S-1,geo.N+1),ylabel=L"Latitude / $\degree$",
			suptitle = L"$p_\omega$" * " ($yr $(monthname(mo)))"
		)
	end

	axs[3].colorbar(c1,label=L"$\mu$ / hPa")
	axs[3].colorbar(c2,label=L"$\Delta$ / hPa")
	
	fig.savefig(projectdir("figures","figS2c-wrfvsera-$(Dates.format(Date(yr,mo),dateformat"yyyymm")).png"),transparent=false,dpi=400)
	load(projectdir("figures","figS2c-wrfvsera-$(Dates.format(Date(yr,mo),dateformat"yyyymm")).png"))
end

# ╔═╡ Cell order:
# ╟─a8431d46-46fe-11ec-2b8d-e39caffdabec
# ╟─b3bc56eb-43a7-4736-bd66-704529911d60
# ╟─9071763e-f6ad-4468-ae48-369307a85263
# ╟─4579f773-e2b4-4b1f-a8c1-485341d451fc
# ╠═7592bde2-b5dc-4b1e-8d51-4a9c33ba25d3
# ╠═c0799dcd-7c88-4f75-ac19-76b615ff4454
# ╠═dc809308-c782-4690-bd79-b4d0d2683cf1
# ╠═538d96e5-ed7b-477a-8f3a-959a9b978ec5
# ╠═9a455f45-592c-43ab-a71c-d1de27bf8552
# ╟─a822f787-214a-4ef3-ba54-b899e890aeff
# ╟─5260b50e-7436-4afa-96b4-3e9ac1d26229
# ╠═9520065b-987c-4b08-8218-3aca80cd4d20
# ╠═696d9f79-c66c-4ff7-9564-aa2fae5cf76b
# ╟─cf67be7e-e020-4a47-8ddc-320f95194c7c
# ╟─d779d4b0-0c3b-4a56-bc44-0ada500c529b
# ╟─19cd09a4-9f58-49cc-ba13-326ce1392413
# ╟─42391f74-64d6-4972-9074-7e923e66c077
# ╟─1b46ee2b-3dfd-460f-b862-000ba7b51617
# ╟─330e270e-4203-4bfc-be20-98f3d9ae89fa
# ╟─b74efe40-6288-4b27-b757-5d0771f2552e

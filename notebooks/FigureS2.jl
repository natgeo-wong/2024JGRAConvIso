### A Pluto.jl notebook ###
# v0.20.5

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
npd = IMERGFinalHH(start=Date(2019,8,1),stop=Date(2020,12,31),path=datadir())

# ╔═╡ c0799dcd-7c88-4f75-ac19-76b615ff4454
e5ds = ERA5Hourly(start=Date(2019,8,1),stop=Date(2020,12,31),path=datadir())

# ╔═╡ 538d96e5-ed7b-477a-8f3a-959a9b978ec5
geo = GeoRegion("OTREC_wrf_d02",path=srcdir())

# ╔═╡ 9a455f45-592c-43ab-a71c-d1de27bf8552
egeo = ERA5Region(geo)

# ╔═╡ a822f787-214a-4ef3-ba54-b899e890aeff
nlsd = getLandSea(npd,geo)

# ╔═╡ 5260b50e-7436-4afa-96b4-3e9ac1d26229
elsd = getLandSea(e5ds,egeo)

# ╔═╡ 0e89ce44-1788-4744-8170-cccb4aecd220
begin
	prcp  = zeros(length(nlsd.lon),length(nlsd.lat))
	dtvec = npd.start : Day(1) : npd.stop
	for idt in dtvec
		ids = read(npd,geo,idt)
		prcp[:,:] += mean(ids["precipitation"][:,:,:],dims=3) * 86400
		close(ids)
	end
	prcp ./= length(dtvec)
end

# ╔═╡ d982d95e-e9d8-4e5b-a706-0589a4eb4df8
begin
	wnds = NCDataset(datadir("wrf","regridded","gpm-RAINNC-20190801_20201231.nc"))
	wprcp = dropdims(mean(wnds["RAINNC"][:,:,1:(end-24)],dims=3),dims=3)*24
	close(wnds)
end

# ╔═╡ 5c9367e4-35a7-4e5a-83e6-ff70a108e340
begin
	eprcp = zeros(length(elsd.lon),length(elsd.lat))
	edtvec = e5ds.start : Month(1) : e5ds.stop
	evar = SingleVariable("tp")
	for idt in edtvec
		ids = read(e5ds,evar,egeo,idt)
		eprcp[:,:] += sum(ids[evar.ID][:,:,:],dims=3) * 1000
		close(ids)
	end
	eprcp ./= length(dtvec)
end

# ╔═╡ 6c91263c-a1e2-41f6-83f1-e69e55ef0c38
begin
	weds = NCDataset(datadir("wrf","regridded","era-RAINNC-20190801_20201231.nc"))
	weprcp = dropdims(mean(weds["RAINNC"][:,:,1:(end-24)],dims=3),dims=3)*24
	close(weds)
end

# ╔═╡ b74efe40-6288-4b27-b757-5d0771f2552e
begin
	asp = (geo.N-geo.S+2)/(geo.E-geo.W+2)
	pplt.close(); fig,axs = pplt.subplots(nrows=2,ncols=3,axwidth=1.5,aspect=asp)
	
	c1 = axs[1].pcolormesh(nlsd.lon,nlsd.lat,prcp',levels=2:2:30,cmap="Blues",extend="both")
	axs[2].pcolormesh(nlsd.lon,nlsd.lat,wprcp',levels=2:2:30,cmap="Blues",extend="both")
	c2 = axs[3].pcolormesh(nlsd.lon,nlsd.lat,(wprcp.-prcp)',levels=vcat(-5:-1,-0.5,0.5,1:5)*4,cmap="drywet",extend="both")
	
	axs[4].pcolormesh(elsd.lon,elsd.lat,eprcp',levels=2:2:30,cmap="Blues",extend="both")
	axs[5].pcolormesh(elsd.lon,elsd.lat,weprcp',levels=2:2:30,cmap="Blues",extend="both")
	axs[6].pcolormesh(elsd.lon,elsd.lat,(weprcp.-eprcp)',levels=vcat(-5:-1,-0.5,0.5,1:5)*4,cmap="drywet",extend="both")

	axs[1].format(ultitle="IMERGv7")
	axs[2].format(ultitle="Regridded WRF")
	axs[3].format(ultitle="WRF - IMERGv7")
	axs[4].format(ultitle="ERA5")
	axs[5].format(ultitle="Regridded WRF")
	axs[6].format(ultitle="WRF - ERA5")

	for ax in axs
		ax.plot(x,y,lw=0.5,c="k")
		ax.format(
			xlim=(geo.W-1,geo.E+1),xlabel=L"Longitude / $\degree$",
			ylim=(geo.S-1,geo.N+1),ylabel=L"Latitude / $\degree$",
			suptitle = "(a) Mean Rainfall Rate (2019 Aug - 2020 Dec)"
		)
	end

	axs[3].colorbar(c1,label=L"$\mu$ / mm day$^{-1}$")
	axs[6].colorbar(c2,label=L"$\Delta$ / mm day$^{-1}$")
	
	fig.savefig(projectdir("figures","figS2-imergvswrfvsera.png"),transparent=false,dpi=400)
	load(projectdir("figures","figS2-imergvswrfvsera.png"))
end

# ╔═╡ Cell order:
# ╟─a8431d46-46fe-11ec-2b8d-e39caffdabec
# ╟─b3bc56eb-43a7-4736-bd66-704529911d60
# ╟─9071763e-f6ad-4468-ae48-369307a85263
# ╟─4579f773-e2b4-4b1f-a8c1-485341d451fc
# ╠═7592bde2-b5dc-4b1e-8d51-4a9c33ba25d3
# ╠═c0799dcd-7c88-4f75-ac19-76b615ff4454
# ╠═538d96e5-ed7b-477a-8f3a-959a9b978ec5
# ╠═9a455f45-592c-43ab-a71c-d1de27bf8552
# ╟─a822f787-214a-4ef3-ba54-b899e890aeff
# ╟─5260b50e-7436-4afa-96b4-3e9ac1d26229
# ╟─0e89ce44-1788-4744-8170-cccb4aecd220
# ╠═d982d95e-e9d8-4e5b-a706-0589a4eb4df8
# ╟─5c9367e4-35a7-4e5a-83e6-ff70a108e340
# ╠═6c91263c-a1e2-41f6-83f1-e69e55ef0c38
# ╠═b74efe40-6288-4b27-b757-5d0771f2552e

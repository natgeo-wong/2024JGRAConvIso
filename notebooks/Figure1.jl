### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ e32a00ee-5f32-47a1-a983-91fb77bc5d18
begin
	using Pkg; Pkg.activate()
	using DrWatson
	md"Using DrWatson to ensure reproducibility between different machines ..."
end

# ╔═╡ bec4e6f2-c2ea-421e-8c57-33e1ef90aa21
begin
	@quickactivate "2024GLConvIso"
	using DelimitedFiles
	using ERA5Reanalysis
	using NCDatasets
	using Statistics
	using PyCall, LaTeXStrings
	using PNGFiles, ImageShow

	pplt = pyimport("proplot")

	include(srcdir("common.jl"))
	
	md"Loading modules for the 2024GLConvIso project..."
end

# ╔═╡ 2e7c33da-f8b5-11ec-08f2-2581af96575f
md"
# Figure 1. Theoretical $p_\omega$ Schematic
"

# ╔═╡ 1cfa1b51-5a64-4945-9e61-82a27900f9de
begin
	coast = readdlm(datadir("coast.cst"),comments=true)
	clon  = coast[:,1]
	clat  = coast[:,2]
	md"Preloading coastline data"
end

# ╔═╡ 75925166-d1e8-450e-bae6-29affd25d635
md"
### A. Defining Datasets, Variables and Regions
"

# ╔═╡ a5487828-a442-4a64-b784-65a7319fc90c
e5ds = ERA5Dummy(path=datadir())

# ╔═╡ cb7c6118-e25b-462a-84a5-751ea0682b52
elsd = getLandSea(e5ds,ERA5Region("TRP"))

# ╔═╡ b68195cb-cf2e-4ce4-9999-1d71eacedf6a
md"
### B. Loading ERA5 and GPM Datasets
"

# ╔═╡ acf26a3f-6f2c-4bfa-a0b4-1d0379cc47fe
begin
	wpds = NCDataset(datadir("p_wwgt-climatology.nc"))
	lon  = wpds["longitude"][:]
	lat  = wpds["latitude"][:]
	wp_yr = wpds["p_wwgt"][:,:] / 100
	sp_yr = wpds["sp"][:,:] / 100
	close(wpds)
	wp_yr[elsd.z.>500] .= NaN
	md"Loading and filtering w-weighted pressure and sigma data and counts ..."
end

# ╔═╡ 7c954c2b-bfdc-4efd-a3b4-8038e66e1ee5
begin
	wds = NCDataset(datadir("fig1wrfdata.nc"))
	μq  = wds["QVAPOR"][:]
	μqh = wds["HDO_QVAPOR"][:]
	μp  = wds["P"][:]
	μδ  = (μqh ./ μq .- 1) * 1000
	close(wds)
	md"Loading the WRF 2D and 3D datasets"
end

# ╔═╡ ee25250d-ddad-4242-9536-d3d0f72e0f85
begin
	cds = NCDataset(datadir("fig1clouddata.nc"))
	xx  = cds["x"][:]
	zz  = cds["z"][:]
	qc  = cds["qc"][:,:]
	close(cds)
end

# ╔═╡ 00f069bb-0834-4859-94fa-7b0c146e3a11
md"
### C. Plotting the climatology
"

# ╔═╡ 4342521a-6a0a-4d30-a5f3-975514d9006f
begin
	pplt.close(); fig,axs = pplt.subplots([
		[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
		[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
		[2,2,2,2,2,2,2,2,2,0,0,0,3,3,3,3,3,3,3,3,3],
		[2,2,2,2,2,2,2,2,2,0,0,0,3,3,3,3,3,3,3,3,3],
		[2,2,2,2,2,2,2,2,2,0,0,0,3,3,3,3,3,3,3,3,3],
		[2,2,2,2,2,2,2,2,2,0,0,0,3,3,3,3,3,3,3,3,3],
		[4,4,4,4,4,4,4,4,4,6,0,7,5,5,5,5,5,5,5,5,5],
		[4,4,4,4,4,4,4,4,4,6,0,7,5,5,5,5,5,5,5,5,5],
		[4,4,4,4,4,4,4,4,4,6,0,7,5,5,5,5,5,5,5,5,5],
	],aspect=7,axwidth=5,wspace=0,sharex=0)
	
	c1 = axs[1].pcolormesh(elsd.lon,elsd.lat,wp_yr',levels=(8:0.5:16).*50,extend="both",cmap="drywet_r")

	for ii in 1 : 3
		axs[ii].pcolormesh(elsd.lon,elsd.lat,wp_yr',levels=(8:0.5:16).*50,extend="both",cmap="drywet_r")
		axs[ii].plot(clon,clat,c="k",lw=0.5)
		axs[ii].format(ylabel=L"Latitude / $\degree$",ylocator=-10:10:20,yminorlocator=-10:5:20,xminorlocator=90:5:300)
	end

	axs[1].plot([100,130,130,100,100],[-5,-5,15,15,-5],c="k",linestyle="--")
	axs[1].plot([260,290,290,260,260],[-5,-5,15,15,-5],c="k",linestyle="--")
	axs[1].text(133,10,"(b)")
	axs[1].text(251,-3,"(c)")

	axs[2].plot([0,360],[0,0],c="r",linestyle="--")
	axs[3].plot([0,360],[5,5],c="r",linestyle="--")
	
	axs[1].format(xlim=(90,300),xlocator=60:30:330,ylim=(-10,20),ltitle="(a) Tropical Pacific Basin")
	axs[2].format(xlim=(100,130),xlocator=100:5:130,ylim=(-5,15),ylocator=-15:5:30,yminorlocator=-15:20,xminorlocator=100:130,ltitle="(b) West Pacific",ultitle="(i)")
	axs[3].format(xlim=(260,290),xlocator=260:5:290,ylim=(-5,15),ylocator=-15:5:20,yminorlocator=-15:20,xminorlocator=260:290,ltitle="(c) East Pacific",ultitle="(i)")

	c2 = axs[4].contourf(elsd.lon,μp/100,repeat(μδ,outer=(1,1440)),levels=-600:50:-150,cmap="viridis",extend="both",cmap_kw=Dict("left"=>0.2))
	axs[4].fill_between(elsd.lon,1000,sp_yr[:,121],c="brown")
	axs[4].plot([0,360],ones(2)*10. ^2.65,c="r",linestyle=":")
	axs[4].format(xlim=(100,130),xlabel=L"Longitude / $\degree$",xlocator=100:5:130,xminorlocator=100:130,ultitle="(ii)")

	axs[5].contourf(elsd.lon,μp/100,repeat(μδ,outer=(1,1440)),levels=-600:50:-150,cmap="viridis",extend="both",cmap_kw=Dict("left"=>0.2))
	axs[5].fill_between(elsd.lon,1000,sp_yr[:,101],c="brown")
	axs[5].plot([0,360],ones(2)*10. ^2.8,c="r",linestyle=":")
	axs[5].format(xlim=(260,290),xlabel=L"Longitude / $\degree$",xlocator=260:5:290,xminorlocator=260:290,ultitle="(ii)")

	axs[6].plot(1.5sin.((0:0.01:1)*pi),10. .^(LinRange(2.3,3,101)),c="k")
	axs[6].plot([0,0],10. .^[2,2.3],c="k")
	axs[6].plot([-5,5],ones(2)*10. ^2.65,c="r",linestyle=":")
	axs[6].format(ultitle="(iii)",xticklabels=[" "])
	
	axs[7].plot(0.4sin.((0:0.01:1)*pi).-sin.((0:0.01:1)*2pi),10. .^(LinRange(2.3,3,101)),c="k")
	axs[7].plot([0,0],10. .^[2,2.3],c="k")
	axs[7].plot([-5,5],ones(2)*10. ^2.8,c="r",linestyle=":")
	axs[7].format(urtitle="(iii)",xticklabels=[" "])

	for ii in 4 : 5
		axs[ii].format(ylim=(1000,100),yscale="log",xlocator=0:5:360,ylabel="Pressure / hPa",xlabel=L"Longitude / $\degree$")
	end

	for ii in 6 : 7
		axs[ii].format(xlim=(-2,2),xticks=[0],yloc="none",xloc="b",xlabel=L"w")
	end

	fig.colorbar(c1,locator=400:50:800,minorlocator=(5:0.5:15).*50,label=L"$p_\omega$ / hPa",length=0.7,rows=(1,6))
	fig.colorbar(c2,locator=[-600,-150],label=L"$\delta_H$",ticklabels=[L"-",L"+"],rows=(7,9))
	fig.savefig(projectdir("figures","fig1-pomegaclimatology.png"),transparent=false,dpi=400)
	load(projectdir("figures","fig1-pomegaclimatology.png"))
end

# ╔═╡ cc7368e3-9591-4bf6-a907-67c585106a8c
md"
### D. Loading some 2D Cloud images
"

# ╔═╡ c07e34e0-01f5-4e0e-afa1-bd3d5f157a87
begin
	pplt.close(); ftmp,atmp = pplt.subplots(ncols=2,aspect=2,axwidth=4,sharex=0)
	
	atmp[1].contourf(xx[8800:9600]/1000,zz/1000,qc[8800:9600,:]',cmap="greys")
	atmp[2].contourf(xx[15950:16250]/1000,zz/1000,qc[15950:16250,:]',cmap="greys")

	for ax in atmp
		ax.format(ylim=(0,12),grid=false,xticks="null",yticks="null")
	end

	ftmp.savefig(projectdir("figures","fig1-clouds.png"),transparent=true,dpi=400)
	load(projectdir("figures","fig1-clouds.png"))
end

# ╔═╡ Cell order:
# ╟─2e7c33da-f8b5-11ec-08f2-2581af96575f
# ╟─e32a00ee-5f32-47a1-a983-91fb77bc5d18
# ╟─bec4e6f2-c2ea-421e-8c57-33e1ef90aa21
# ╟─1cfa1b51-5a64-4945-9e61-82a27900f9de
# ╟─75925166-d1e8-450e-bae6-29affd25d635
# ╟─a5487828-a442-4a64-b784-65a7319fc90c
# ╟─cb7c6118-e25b-462a-84a5-751ea0682b52
# ╟─b68195cb-cf2e-4ce4-9999-1d71eacedf6a
# ╠═acf26a3f-6f2c-4bfa-a0b4-1d0379cc47fe
# ╟─7c954c2b-bfdc-4efd-a3b4-8038e66e1ee5
# ╟─ee25250d-ddad-4242-9536-d3d0f72e0f85
# ╟─00f069bb-0834-4859-94fa-7b0c146e3a11
# ╟─4342521a-6a0a-4d30-a5f3-975514d9006f
# ╟─cc7368e3-9591-4bf6-a907-67c585106a8c
# ╟─c07e34e0-01f5-4e0e-afa1-bd3d5f157a87

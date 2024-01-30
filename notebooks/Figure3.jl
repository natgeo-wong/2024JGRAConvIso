### A Pluto.jl notebook ###
# v0.19.37

using Markdown
using InteractiveUtils

# ╔═╡ 9802aaa7-3f9a-47b7-b6ab-90c4f39b7335
begin
	using Pkg; Pkg.activate()
	using DrWatson
	md"Using DrWatson to ensure reproducibility between different machines ..."
end

# ╔═╡ b62ba51a-b7a8-433c-84dd-bd7d221ffa3c
begin
	@quickactivate "2024GLConvIso"
	using Dates
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

# ╔═╡ fc7b6caa-6ced-11ec-0701-6f55729e22dc
md"
# Figure 3. Binning $\delta$ against $p_\omega$ and $P$
"

# ╔═╡ a7e28ecc-6808-49d8-838d-5585336c01aa
md"
### A. ERA5 Datasets
"

# ╔═╡ 3d568d03-50b2-44ca-ae3c-c8c383bf3e9f
e5dy = ERA5Daily(start=Date(2013),stop=Date(2021),path=datadir())

# ╔═╡ afb4338b-4f4d-406f-af61-630ed6f85297
evar = SingleVariable("p_wwgt")

# ╔═╡ 21894078-46cc-43c8-a465-047d5318187c
ereg = ERA5Region(GeoRegion("OTREC"))

# ╔═╡ 65a509df-0943-4054-ae0a-216e405aad48
lsd = getLandSea(e5dy,ereg)

# ╔═╡ a1007189-6c4d-4418-aefa-078c0090def5
md"
### B. Plotting functions and the like
"

# ╔═╡ 2c696d1a-efd3-4b41-bc25-d0c0e8ac2543
function threshold!(
	bin,num,prc;
	numthresh = 5
)

	bin[num.<numthresh] .= 0
	prc[num.<numthresh] .= 0
	# num[num.<numthresh] .= 0

	return
	
end

# ╔═╡ 56109fca-661d-4594-87b4-677511437a2e
function plotbin!(
	axesnum,ii,
	rainbin,wgtpbin,
	iibin, iiprc, iinum,
	lvls;
	returncinfo = true
)

	c1 = axesnum[2*ii-1].pcolormesh(
		rainbin,wgtpbin,(iibin./iiprc)',
		cmap="viridis",levels=lvls,extend="both"
	)
	c2 = axesnum[2*ii].pcolormesh(
		rainbin,wgtpbin,iinum',
		cmap="fire",levels=0:10:100,extend="both"
	)

	ix = axesnum[2*ii-1].panel("l",width="0.6em",space=0)
	ix.pcolormesh(
		[0,1],wgtpbin,(sum(iibin,dims=1)./sum(iiprc,dims=1))',
		cmap="viridis",levels=lvls,extend="both"
	)
	ix.format(xlocator=[])

	ix = axesnum[2*ii-1].panel("t",width="0.6em",space=0)
	ix.pcolormesh(
		rainbin,[0,1],(sum(iibin,dims=2)./sum(iiprc,dims=2))',
		cmap="viridis",levels=lvls,extend="both"
	)
	ix.format(ylocator=[],xlabel=L"$P$ / kg m$^{-2}$ day$^{-1}$")

	ix = axesnum[2*ii].panel("r",width="0.6em",space=0)
	ix.pcolormesh(
		[0,1],wgtpbin,sum(iinum,dims=1)',
		cmap="fire",levels=0:10:100,extend="both"
	)
	ix.format(xlocator=[])

	ix = axesnum[2*ii].panel("t",width="0.6em",space=0)
	ix.pcolormesh(
		rainbin,[0,1],sum(iinum,dims=2)',
		cmap="fire",levels=0:10:100,extend="both"
	)
	ix.format(ylocator=[],xlabel=L"$P$ / kg m$^{-2}$ day$^{-1}$")

	if returncinfo
		return c1,c2
	else
		return
	end

end

# ╔═╡ b266abf5-d287-4ca0-945e-2c771529fd41
function axesformat!(axesnum)

	for ax in axesnum
		ax.format(
			ylabel=L"$p_\omega$ / hPa",
			xlabel=L"Rainfall Rate / mm day$^{-1}$",
			ylim=(950,100),xlim=(0,100),ylocator=200:200:800,xlocator=0:50:100
		)
	end
	
	axesnum[2].format(urtitle="(a) All Stations")
	axesnum[4].format(urtitle="(b) Colombia")
	axesnum[6].format(urtitle="(c) Costa Rica")
	axesnum[8].format(urtitle="(d) San Andres")
	axesnum[10].format(urtitle="(e) Buenaventura")
	axesnum[10].text(-23,300,"Bahia Solano",fontsize=10)
	axesnum[12].format(urtitle="(f) Quibdo")

	return

end

# ╔═╡ 8927283a-c4e2-4c3a-b32f-850391eef9c7
md"
### C. Loading and Comparing Daily Station Data
"

# ╔═╡ f4597370-4bc8-4763-8030-3a4dc89533b6
begin
	infody  = stninfody(); ndy = size(infody,1)
	md"Loading daily station location information ..."
end

# ╔═╡ 2c4f6363-7b5a-40a0-aa47-251c53ae7367
begin
	ds = NCDataset(datadir("processed.nc"))
	time  = ds["time"][:]; it = time .> Date(2019,1,31)
	time  = time[it]
	prcps = ds["prcps"][it,:]
	δ18Oμ = ds["δ18Oμ"][it,:]
	δHDOμ = ds["δ2Hμ"][it,:]
	nstn = size(prcps,2)
	close(ds)
end

# ╔═╡ b2dc480b-b6fc-471a-9724-1034415f4cfc
begin
	rbin = 0 : 10 : 100;    rpnt = (rbin[1:(end-1)] .+ rbin[2:end]) / 2
	pbin = 100 : 50 : 950; ppnt = (pbin[1:(end-1)] .+ pbin[2:end]) / 2
	abin = zeros(length(rpnt),length(ppnt)); anum = zeros(length(rpnt),length(ppnt))
	bbin = zeros(length(rpnt),length(ppnt)); bnum = zeros(length(rpnt),length(ppnt))
	cbin = zeros(length(rpnt),length(ppnt)); cnum = zeros(length(rpnt),length(ppnt))
	dbin = zeros(length(rpnt),length(ppnt)); dnum = zeros(length(rpnt),length(ppnt))
	ebin = zeros(length(rpnt),length(ppnt)); enum = zeros(length(rpnt),length(ppnt))
	fbin = zeros(length(rpnt),length(ppnt)); fnum = zeros(length(rpnt),length(ppnt))
	aprc = zeros(length(rpnt),length(ppnt))
	bprc = zeros(length(rpnt),length(ppnt))
	cprc = zeros(length(rpnt),length(ppnt))
	dprc = zeros(length(rpnt),length(ppnt))
	eprc = zeros(length(rpnt),length(ppnt))
	fprc = zeros(length(rpnt),length(ppnt))
	md"Preallocation of arrays ..."
end

# ╔═╡ 3d3859f7-7c57-4bf4-b3d8-0125c10660cb
function binning!(
	binstn,numstn,prcstn,
	rpnt, ppnt;
	ID, days=0, iso = "HDO", bnum = 1
)

	iso = "$(iso)_"
	stnstr = @sprintf("%02d",ID)

	
	geoname = "OTREC_STN$stnstr"
	pwgt,prcp,hvyp = extract(geoname,iso,days)
	nt = length(prcp)

	for it = 1 : nt
		if (prcp[it]>2.5) && !isnan(pwgt[it])
			rind = argmin(abs.(prcp[it].-rpnt))
			pind = argmin(abs.(pwgt[it].-ppnt))
			numstn[rind,pind] += 1
			binstn[rind,pind] += hvyp[it]
			prcstn[rind,pind] += prcp[it]
		end
	end

	for idy = 4 : (length(time)-3)

		dt  = time[idy]; dayii = day(dt)
		ids = read(e5dy,evar,ereg,dt,smooth=true,smoothtime=7,quiet=true)

		prcps_ii = prcps[idy.+(-3:3),:]
		δ18Oμ_ii = δ18Oμ[idy.+(-3:3),:]

		lon = infody[ID,2]; ilon = argmin(abs.(lsd.lon.-lon))
		lat = infody[ID,3]; ilat = argmin(abs.(lsd.lat.-lat))
		prcpsii = prcps_ii[:,ID]; iNaN_prcp = .!isnan.(prcpsii)
		δ18Oμii = δ18Oμ_ii[:,ID]; iNaN_δ18O = .!isnan.(δ18Oμii)
		iNaN = iNaN_prcp .& iNaN_δ18O
		if sum(iNaN) > thresholdNaN
			prcpμii = mean(prcpsii[iNaN])
			prcptot = sum(prcpsii[iNaN])
			δ18Oμii = sum(δ18Oμii[iNaN] .* prcpsii[iNaN])
			pwwgtii = ids["p_wwgt"][ilon,ilat,dayii] / 100
			if ismissing(pwwgtii); pwwgtii = NaN end
	
			if !isnan(prcpμii) && !isnan(pwwgtii)
				rind = argmin(abs.(prcpμii.-rpnt))
				pind = argmin(abs.(pwwgtii.-ppnt))
				abin[rind,pind] += δ18Oμii; anum[rind,pind] += 1
				bbin[rind,pind] += δ18Oμii; bnum[rind,pind] += 1
				aprc[rind,pind] += prcptot; bprc[rind,pind] += prcptot
			end
		end

	end

	return

end

# ╔═╡ 492b98ea-e5d8-40e3-ae89-df68d183c100
md"Now, we do a 7-day moving average ..."

# ╔═╡ 189bb0b7-d3ae-4aa7-b9e3-1ee5dcc63e4d
ithr_7dy = 4

# ╔═╡ bbcf49bb-4be0-463d-87ea-83f674a5729d
begin
	pplt.close()
	f3,a3 = pplt.subplots(ncols=6,nrows=2,aspect=1/2,axwidth=0.75,wspace=[0,2,0,2,0])

	lvls = -10 : -1
	abin .= 0; bbin .= 0; cbin .= 0; dbin .= 0; ebin .= 0; fbin .= 0
	anum .= 0; bnum .= 0; cnum .= 0; dnum .= 0; enum .= 0; fnum .= 0
	aprc .= 0; bprc .= 0; cprc .= 0; dprc .= 0; eprc .= 0; fprc .= 0

	for idy = 4 : (length(time)-3)

		dt  = time[idy]; dayii = day(dt)
		ids = read(e5dy,evar,ereg,dt,smooth=true,smoothtime=7,quiet=true)

		prcps_ii = prcps[idy.+(-3:3),:]
		δ18Oμ_ii = δ18Oμ[idy.+(-3:3),:]

		for istn in [1]

			lon = infody[istn,2]; ilon = argmin(abs.(lsd.lon.-lon))
			lat = infody[istn,3]; ilat = argmin(abs.(lsd.lat.-lat))
			prcpsii = prcps_ii[:,istn]; iNaN_prcp = .!isnan.(prcpsii)
			δ18Oμii = δ18Oμ_ii[:,istn]; iNaN_δ18O = .!isnan.(δ18Oμii)
			iNaN = iNaN_prcp .& iNaN_δ18O
			if sum(iNaN) > ithr_7dy
				prcpμii = mean(prcpsii[iNaN])
				prcptot = sum(prcpsii[iNaN])
				δ18Oμii = sum(δ18Oμii[iNaN] .* prcpsii[iNaN])
				pwwgtii = ids["p_wwgt"][ilon,ilat,dayii] / 100
				if ismissing(pwwgtii); pwwgtii = NaN end
	
				if !isnan(prcpμii) && !isnan(pwwgtii)
					rind = argmin(abs.(prcpμii.-rpnt))
					pind = argmin(abs.(pwwgtii.-ppnt))
					abin[rind,pind] += δ18Oμii; anum[rind,pind] += 1
					bbin[rind,pind] += δ18Oμii; bnum[rind,pind] += 1
					aprc[rind,pind] += prcptot; bprc[rind,pind] += prcptot
				end
			end

		end

		for istn in [3,4]

			lon = infody[istn,2]; ilon = argmin(abs.(lsd.lon.-lon))
			lat = infody[istn,3]; ilat = argmin(abs.(lsd.lat.-lat))
			prcpsii = prcps_ii[:,istn]; iNaN_prcp = .!isnan.(prcpsii)
			δ18Oμii = δ18Oμ_ii[:,istn]; iNaN_δ18O = .!isnan.(δ18Oμii)
			iNaN = iNaN_prcp .& iNaN_δ18O
			if sum(iNaN) > ithr_7dy
				prcpμii = mean(prcpsii[iNaN])
				prcptot = sum(prcpsii[iNaN])
				δ18Oμii = sum(δ18Oμii[iNaN] .* prcpsii[iNaN])
				pwwgtii = ids["p_wwgt"][ilon,ilat,dayii] / 100
				if ismissing(pwwgtii); pwwgtii = NaN end
	
				if !isnan(prcpμii) && !isnan(pwwgtii)
					rind = argmin(abs.(prcpμii.-rpnt))
					pind = argmin(abs.(pwwgtii.-ppnt))
					abin[rind,pind] += δ18Oμii; anum[rind,pind] += 1
					cbin[rind,pind] += δ18Oμii; cnum[rind,pind] += 1
					aprc[rind,pind] += prcptot; cprc[rind,pind] += prcptot
				end
			end

		end

		for istn in [2]

			lon = infody[istn,2]; ilon = argmin(abs.(lsd.lon.-lon))
			lat = infody[istn,3]; ilat = argmin(abs.(lsd.lat.-lat))
			prcpsii = prcps_ii[:,istn]; iNaN_prcp = .!isnan.(prcpsii)
			δ18Oμii = δ18Oμ_ii[:,istn]; iNaN_δ18O = .!isnan.(δ18Oμii)
			iNaN = iNaN_prcp .& iNaN_δ18O
			if sum(iNaN) > ithr_7dy
				prcpμii = mean(prcpsii[iNaN])
				prcptot = sum(prcpsii[iNaN])
				δ18Oμii = sum(δ18Oμii[iNaN] .* prcpsii[iNaN])
				pwwgtii = ids["p_wwgt"][ilon,ilat,dayii] / 100
				if ismissing(pwwgtii); pwwgtii = NaN end
	
				if !isnan(prcpμii) && !isnan(pwwgtii)
					rind = argmin(abs.(prcpμii.-rpnt))
					pind = argmin(abs.(pwwgtii.-ppnt))
					abin[rind,pind] += δ18Oμii; anum[rind,pind] += 1
					dbin[rind,pind] += δ18Oμii; dnum[rind,pind] += 1
					aprc[rind,pind] += prcptot; dprc[rind,pind] += prcptot
				end
			end

		end

		for istn in 5 : 12

			lon = infody[istn,2]; ilon = argmin(abs.(lsd.lon.-lon))
			lat = infody[istn,3]; ilat = argmin(abs.(lsd.lat.-lat))
			prcpsii = prcps_ii[:,istn]; iNaN_prcp = .!isnan.(prcpsii)
			δ18Oμii = δ18Oμ_ii[:,istn]; iNaN_δ18O = .!isnan.(δ18Oμii)
			iNaN = iNaN_prcp .& iNaN_δ18O
			if sum(iNaN) > ithr_7dy
				prcpμii = mean(prcpsii[iNaN])
				prcptot = sum(prcpsii[iNaN])
				δ18Oμii = sum(δ18Oμii[iNaN] .* prcpsii[iNaN])
				pwwgtii = ids["p_wwgt"][ilon,ilat,dayii] / 100
				if ismissing(pwwgtii); pwwgtii = NaN end
	
				if !isnan(prcpμii) && !isnan(pwwgtii)
					rind = argmin(abs.(prcpμii.-rpnt))
					pind = argmin(abs.(pwwgtii.-ppnt))
					abin[rind,pind] += δ18Oμii; anum[rind,pind] += 1
					ebin[rind,pind] += δ18Oμii; enum[rind,pind] += 1
					aprc[rind,pind] += prcptot; eprc[rind,pind] += prcptot
				end
			end

		end

		for istn in 1 : 4

			lon = infody[istn,2]; ilon = argmin(abs.(lsd.lon.-lon))
			lat = infody[istn,3]; ilat = argmin(abs.(lsd.lat.-lat))
			prcpsii = prcps_ii[:,istn]; iNaN_prcp = .!isnan.(prcpsii)
			δ18Oμii = δ18Oμ_ii[:,istn]; iNaN_δ18O = .!isnan.(δ18Oμii)
			iNaN = iNaN_prcp .& iNaN_δ18O
			if sum(iNaN) > ithr_7dy
				prcpμii = mean(prcpsii[iNaN])
				prcptot = sum(prcpsii[iNaN])
				δ18Oμii = sum(δ18Oμii[iNaN] .* prcpsii[iNaN])
				pwwgtii = ids["p_wwgt"][ilon,ilat,dayii] / 100
				if ismissing(pwwgtii); pwwgtii = NaN end
	
				if !isnan(prcpμii) && !isnan(pwwgtii)
					rind = argmin(abs.(prcpμii.-rpnt))
					pind = argmin(abs.(pwwgtii.-ppnt))
					fbin[rind,pind] += δ18Oμii; fnum[rind,pind] += 1
					fprc[rind,pind] += prcptot
				end
			end

		end

		close(ids)

	end

	c3_1,c3_2 = 
	plotbin!(a3,1,rbin,pbin,abin,aprc,anum,-10:0,returncinfo=true)
	plotbin!(a3,2,rbin,pbin,fbin,fprc,fnum,-10:0,returncinfo=true)
	plotbin!(a3,3,rbin,pbin,ebin,eprc,enum,-10:0,returncinfo=true)
	plotbin!(a3,4,rbin,pbin,bbin,bprc,bnum,-10:0,returncinfo=true)
	plotbin!(a3,5,rbin,pbin,cbin,cprc,cnum,-10:0,returncinfo=true)
	plotbin!(a3,6,rbin,pbin,dbin,dprc,dnum,-10:0,returncinfo=true)

	axesformat!(a3)
	a3[1].format(suptitle="7-Day WRF Moving Average")
	f3.colorbar(c3_1,label=L"$\delta^{18}$O / $\perthousand$",locator=-70:2:0,minorlocator=-90:0,row=[1])
	f3.colorbar(c3_2,label="Number of Observations",row=[2])

	f3.savefig(projectdir("figures","fig3-observeddepletion.png"),transparent=false,dpi=400)
	load(projectdir("figures","fig3-observeddepletion.png"))
end

# ╔═╡ Cell order:
# ╟─fc7b6caa-6ced-11ec-0701-6f55729e22dc
# ╟─9802aaa7-3f9a-47b7-b6ab-90c4f39b7335
# ╟─b62ba51a-b7a8-433c-84dd-bd7d221ffa3c
# ╟─a7e28ecc-6808-49d8-838d-5585336c01aa
# ╟─3d568d03-50b2-44ca-ae3c-c8c383bf3e9f
# ╟─afb4338b-4f4d-406f-af61-630ed6f85297
# ╟─21894078-46cc-43c8-a465-047d5318187c
# ╟─65a509df-0943-4054-ae0a-216e405aad48
# ╟─a1007189-6c4d-4418-aefa-078c0090def5
# ╟─3d3859f7-7c57-4bf4-b3d8-0125c10660cb
# ╟─2c696d1a-efd3-4b41-bc25-d0c0e8ac2543
# ╟─56109fca-661d-4594-87b4-677511437a2e
# ╟─b266abf5-d287-4ca0-945e-2c771529fd41
# ╟─8927283a-c4e2-4c3a-b32f-850391eef9c7
# ╟─f4597370-4bc8-4763-8030-3a4dc89533b6
# ╟─2c4f6363-7b5a-40a0-aa47-251c53ae7367
# ╟─b2dc480b-b6fc-471a-9724-1034415f4cfc
# ╟─492b98ea-e5d8-40e3-ae89-df68d183c100
# ╠═189bb0b7-d3ae-4aa7-b9e3-1ee5dcc63e4d
# ╠═bbcf49bb-4be0-463d-87ea-83f674a5729d

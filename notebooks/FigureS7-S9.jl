### A Pluto.jl notebook ###
# v0.19.38

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

# ╔═╡ 2ce41596-200d-4b42-8d54-328f495e81fb
function calculatebufferweights(shiftsteps)

    buffer = Int(ceil((shiftsteps-1)/2))
    weights = ones(buffer*2+1)
    if buffer >= (shiftsteps/2)
        weights[1] = 0.5
        weights[end] = 0.5
    end
    weights /= shiftsteps
    return buffer,weights

end

# ╔═╡ ce7187b0-9683-4388-888b-beeec3888cc2
function smooth(data::AbstractVector,days,NaN_threshold;dosum = false)

	buffer,weights = calculatebufferweights(days)

	ndt = length(data)
	ndata = fill(NaN,ndt)
	smth  = zeros(1+buffer*2)

	for ii in (1+buffer) : (ndt-buffer)

		for ismth = 0 : (buffer*2)
			smth[ismth+1] = data[ii+ismth-buffer] * weights[ismth+1]
		end
		if sum(.!isnan.(smth)) > NaN_threshold
			iinan = .!isnan.(smth)
			ndata[ii] = sum(smth[iinan]) / sum(weights[iinan])
			if dosum; ndata[ii] *= sum(iinan)  end
		end
		

	end

	return ndata

end

# ╔═╡ 3d3859f7-7c57-4bf4-b3d8-0125c10660cb
function binning!(
	binstn,numstn,prcstn,
	rpnt, ppnt;
	ID, days=0, iso = "2H", NaN_threshold = days/3*2
)

	ds = NCDataset(datadir("processed.nc"))
	time  = ds["time"][:]; it = time .> Date(2019,1,31)
	time  = time[it]
	prcps = ds["prcps"][it,ID]
	δhvyμ = ds["δ$(iso)μ"][it,ID]
	nstn = size(prcps,2)
	close(ds)

	prcps_ii = smooth(prcps,days,NaN_threshold)
	prcpt_ii = smooth(prcps,days,NaN_threshold,dosum=true)
	δhvyμ_ii = smooth(δhvyμ .* prcps,days,NaN_threshold,dosum=true)

	ids = NCDataset(datadir("processed-p_wwgt.nc"))
	pwwgt = nomissing(ids["pω$(@sprintf("%02d",days))"][:,ID],NaN) / 100
	close(ids)

	for idy = 1 : length(time)

		dt  = time[idy]; dayii = day(dt)

		prcpμii = prcps_ii[idy]
		prcptii = prcpt_ii[idy]
		δhvyμii = δhvyμ_ii[idy]
		pwwgtii = pwwgt[idy]
	
		if !isnan(prcpμii) && !isnan(δhvyμii) && !isnan(pwwgtii)
			rind = argmin(abs.(prcpμii.-rpnt))
			pind = argmin(abs.(pwwgtii.-ppnt))
			binstn[rind,pind] += δhvyμii
			numstn[rind,pind] += 1
			prcstn[rind,pind] += prcptii
		end

	end

	return

end

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
	returncinfo = true,
	doalpha = false
)

	tmpbin = sum(iibin,dims=1)
	tmpprc = sum(iiprc,dims=1)
	tmpnum = sum(iinum,dims=1)
	threshold!(tmpbin,tmpnum,tmpprc)
	ix = axesnum[2*ii-1].panel("l",width="0.6em",space=0)
	ix.pcolormesh(
		[0,1],wgtpbin,(tmpbin./tmpprc)',
		cmap="viridis",levels=lvls,extend="both"
	)
	ix.format(xlocator=[])

	ix = axesnum[2*ii].panel("r",width="0.6em",space=0)
	ix.pcolormesh(
		[0,1],wgtpbin,tmpnum',
		cmap="fire",levels=0:5:50,extend="both"
	)
	ix.format(xlocator=[])

	tmpbin = sum(iibin,dims=2)
	tmpprc = sum(iiprc,dims=2)
	tmpnum = sum(iinum,dims=2)
	threshold!(tmpbin,tmpnum,tmpprc)
	ix = axesnum[2*ii-1].panel("t",width="0.6em",space=0)
	ix.pcolormesh(
		rainbin,[0,1],(tmpbin./tmpprc)',
		cmap="viridis",levels=lvls,extend="both"
	)
	ix.format(ylocator=[])

	ix = axesnum[2*ii].panel("t",width="0.6em",space=0)
	ix.pcolormesh(
		rainbin,[0,1],tmpnum',
		cmap="fire",levels=0:5:50,extend="both"
	)
	ix.format(ylocator=[])

	threshold!(iibin,iinum,iiprc)
	c1 = axesnum[2*ii-1].pcolormesh(
		rainbin,wgtpbin,(iibin./iiprc)',
		cmap="viridis",levels=lvls,extend="both"
	)
	c2 = axesnum[2*ii].pcolormesh(
		rainbin,wgtpbin,iinum',
		cmap="fire",levels=0:5:50,extend="both"
	)

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

# ╔═╡ b2dc480b-b6fc-471a-9724-1034415f4cfc
begin
	rbin = 0 : 10 : 100;   rpnt = (rbin[1:(end-1)] .+ rbin[2:end]) / 2
	pbin = 100 : 50 : 950; ppnt = (pbin[1:(end-1)] .+ pbin[2:end]) / 2
	nr = length(rpnt); np = length(ppnt)
	abin = zeros(nr,np); anum = zeros(nr,np); aprc = zeros(nr,np)
	bbin = zeros(nr,np); bnum = zeros(nr,np); bprc = zeros(nr,np)
	cbin = zeros(nr,np); cnum = zeros(nr,np); cprc = zeros(nr,np)
	dbin = zeros(nr,np); dnum = zeros(nr,np); dprc = zeros(nr,np)
	ebin = zeros(nr,np); enum = zeros(nr,np); eprc = zeros(nr,np)
	fbin = zeros(nr,np); fnum = zeros(nr,np); fprc = zeros(nr,np)
	md"Preallocation of arrays ..."
end

# ╔═╡ bbcf49bb-4be0-463d-87ea-83f674a5729d
begin
	pplt.close()
	f1,a1 = pplt.subplots(ncols=6,nrows=2,aspect=1/2,axwidth=0.75,wspace=[0,2,0,2,0])

	lvls = -10 : -1

	abin .= 0; aprc .= 0; anum .= 0
	bbin .= 0; bprc .= 0; bnum .= 0
	cbin .= 0; cprc .= 0; cnum .= 0
	dbin .= 0; dprc .= 0; dnum .= 0
	ebin .= 0; eprc .= 0; enum .= 0
	fbin .= 0; fprc .= 0; fnum .= 0

	for istn = 1 : 12
		binning!(abin,anum,aprc,rpnt,ppnt,ID=istn,days=7)
	end
	for istn = 1 : 4
		binning!(bbin,bnum,bprc,rpnt,ppnt,ID=istn,days=7)
	end
	for istn = 5 : 12
		binning!(cbin,cnum,cprc,rpnt,ppnt,ID=istn,days=7)
	end
	for istn = 1
		binning!(dbin,dnum,dprc,rpnt,ppnt,ID=istn,days=7)
	end
	for istn = 3 : 4
		binning!(ebin,enum,eprc,rpnt,ppnt,ID=istn,days=7)
	end
	for istn = 2
		binning!(fbin,fnum,fprc,rpnt,ppnt,ID=istn,days=7)
	end

	c1_1,c1_2 = 
	plotbin!(a1,1,rbin,pbin,abin,aprc,anum,-70:5:-10,returncinfo=true)
	plotbin!(a1,2,rbin,pbin,bbin,bprc,bnum,-70:5:-10,returncinfo=true)
	plotbin!(a1,3,rbin,pbin,cbin,cprc,cnum,-70:5:-10,returncinfo=true)
	plotbin!(a1,4,rbin,pbin,dbin,dprc,dnum,-70:5:-10,returncinfo=true)
	plotbin!(a1,5,rbin,pbin,ebin,eprc,enum,-70:5:-10,returncinfo=true)
	plotbin!(a1,6,rbin,pbin,fbin,fprc,fnum,-70:5:-10,returncinfo=true)

	axesformat!(a1)
	a1[1].format(suptitle="7-Day Station Moving Average")
	f1.colorbar(c1_1,label=L"$\delta^{2}$H / $\perthousand$",locator=-70:10:-10,minorlocator=-90:5:0,row=[1])
	f1.colorbar(c1_2,label="Number of Observations",row=[2])

	f1.savefig(projectdir("figures","figS7-observeddepletion-07dayHDO.png"),transparent=false,dpi=400)
	load(projectdir("figures","figS7-observeddepletion-07dayHDO.png"))
end

# ╔═╡ d7ee8f37-49e3-4878-92e2-c70cc773ab28
begin
	pplt.close()
	f2,a2 = pplt.subplots(ncols=6,nrows=2,aspect=1/2,axwidth=0.75,wspace=[0,2,0,2,0])

	abin .= 0; aprc .= 0; anum .= 0
	bbin .= 0; bprc .= 0; bnum .= 0
	cbin .= 0; cprc .= 0; cnum .= 0
	dbin .= 0; dprc .= 0; dnum .= 0
	ebin .= 0; eprc .= 0; enum .= 0
	fbin .= 0; fprc .= 0; fnum .= 0

	for istn = 1 : 12
		binning!(abin,anum,aprc,rpnt,ppnt,ID=istn,iso="18O",days=30)
	end
	for istn = 1 : 4
		binning!(bbin,bnum,bprc,rpnt,ppnt,ID=istn,iso="18O",days=30)
	end
	for istn = 5 : 12
		binning!(cbin,cnum,cprc,rpnt,ppnt,ID=istn,iso="18O",days=30)
	end
	for istn = 1
		binning!(dbin,dnum,dprc,rpnt,ppnt,ID=istn,iso="18O",days=30)
	end
	for istn = 2
		binning!(fbin,fnum,fprc,rpnt,ppnt,ID=istn,iso="18O",days=30)
	end
	for istn = 3 : 4
		binning!(ebin,enum,eprc,rpnt,ppnt,ID=istn,iso="18O",days=30)
	end

	c2_1,c2_2 = 
	plotbin!(a2,1,rbin,pbin,abin,aprc,anum,-10:-2,returncinfo=true)
	plotbin!(a2,2,rbin,pbin,bbin,bprc,bnum,-10:-2,returncinfo=true)
	plotbin!(a2,3,rbin,pbin,cbin,cprc,cnum,-10:-2,returncinfo=true)
	plotbin!(a2,4,rbin,pbin,dbin,dprc,dnum,-10:-2,returncinfo=true)
	plotbin!(a2,5,rbin,pbin,ebin,eprc,enum,-10:-2,returncinfo=true)
	plotbin!(a2,6,rbin,pbin,fbin,fprc,fnum,-10:-2,returncinfo=true)

	axesformat!(a2)
	a2[1].format(suptitle="30-Day Station Moving Average")
	f2.colorbar(c2_1,label=L"$\delta^{18}$O / $\perthousand$",row=[1])
	f2.colorbar(c2_2,label="Number of Observations",row=[2])

	f2.savefig(projectdir("figures","figS8-observeddepletion-30dayO18.png"),transparent=false,dpi=400)
	load(projectdir("figures","figS8-observeddepletion-30dayO18.png"))
end

# ╔═╡ 2b230a91-40a6-401a-92b8-8aadc3f79679
begin
	pplt.close()
	f3,a3 = pplt.subplots(ncols=6,nrows=2,aspect=1/2,axwidth=0.75,wspace=[0,2,0,2,0])

	abin .= 0; aprc .= 0; anum .= 0
	bbin .= 0; bprc .= 0; bnum .= 0
	cbin .= 0; cprc .= 0; cnum .= 0
	dbin .= 0; dprc .= 0; dnum .= 0
	ebin .= 0; eprc .= 0; enum .= 0
	fbin .= 0; fprc .= 0; fnum .= 0

	for istn = 1 : 12
		binning!(abin,anum,aprc,rpnt,ppnt,ID=istn,days=30)
	end
	for istn = 1 : 4
		binning!(bbin,bnum,bprc,rpnt,ppnt,ID=istn,days=30)
	end
	for istn = 5 : 12
		binning!(cbin,cnum,cprc,rpnt,ppnt,ID=istn,days=30)
	end
	for istn = 1
		binning!(dbin,dnum,dprc,rpnt,ppnt,ID=istn,days=30)
	end
	for istn = 2
		binning!(fbin,fnum,fprc,rpnt,ppnt,ID=istn,days=30)
	end
	for istn = 3 : 4
		binning!(ebin,enum,eprc,rpnt,ppnt,ID=istn,days=30)
	end

	c3_1,c3_2 = 
	plotbin!(a3,1,rbin,pbin,abin,aprc,anum,-70:5:-10,returncinfo=true)
	plotbin!(a3,2,rbin,pbin,bbin,bprc,bnum,-70:5:-10,returncinfo=true)
	plotbin!(a3,3,rbin,pbin,cbin,cprc,cnum,-70:5:-10,returncinfo=true)
	plotbin!(a3,4,rbin,pbin,dbin,dprc,dnum,-70:5:-10,returncinfo=true)
	plotbin!(a3,5,rbin,pbin,ebin,eprc,enum,-70:5:-10,returncinfo=true)
	plotbin!(a3,6,rbin,pbin,fbin,fprc,fnum,-70:5:-10,returncinfo=true)

	axesformat!(a3)
	a3[1].format(suptitle="30-Day Station Moving Average")
	f3.colorbar(c3_1,label=L"$\delta^{2}$H / $\perthousand$",locator=-70:10:-10,minorlocator=-90:5:0,row=[1])
	f3.colorbar(c3_2,label="Number of Observations",row=[2])

	f3.savefig(projectdir("figures","figS9-observeddepletion-30dayHDO.png"),transparent=false,dpi=400)
	load(projectdir("figures","figS9-observeddepletion-30dayHDO.png"))
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
# ╟─2ce41596-200d-4b42-8d54-328f495e81fb
# ╟─ce7187b0-9683-4388-888b-beeec3888cc2
# ╟─3d3859f7-7c57-4bf4-b3d8-0125c10660cb
# ╟─2c696d1a-efd3-4b41-bc25-d0c0e8ac2543
# ╟─56109fca-661d-4594-87b4-677511437a2e
# ╟─b266abf5-d287-4ca0-945e-2c771529fd41
# ╟─8927283a-c4e2-4c3a-b32f-850391eef9c7
# ╟─b2dc480b-b6fc-471a-9724-1034415f4cfc
# ╟─bbcf49bb-4be0-463d-87ea-83f674a5729d
# ╟─d7ee8f37-49e3-4878-92e2-c70cc773ab28
# ╟─2b230a91-40a6-401a-92b8-8aadc3f79679

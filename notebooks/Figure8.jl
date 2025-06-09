### A Pluto.jl notebook ###
# v0.20.10

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
	using GeoRegions
	using NCDatasets
	using Printf
	using Statistics
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("ultraplot")

	include(srcdir("common.jl"))
	
	md"Loading modules for the ConvectionIsotopes project..."
end

# ╔═╡ 2e7c33da-f8b5-11ec-08f2-2581af96575f
md"
# 08. A New Index: $p_{q\omega,bl}$
"

# ╔═╡ 59c930cd-5b7f-4047-8660-615148d1bd9f
begin
	infody = stninfody()[:,:]; nstn = size(infody,1)
	md"Loading station location information ..."
end

# ╔═╡ 441f47a7-5757-4b24-8b52-a2877e0f0287
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

# ╔═╡ f1720645-69a8-4f45-a6c1-8c06279d3590
function smooth(data::AbstractVector,days)

	buffer,weights = calculatebufferweights(days)

	ndt = length(data)
	ndata = fill(NaN,ndt)
	smth  = zeros(1+buffer*2)

	for ii in (1+buffer) : (ndt-buffer)

		for ismth = 0 : (buffer*2)
			smth[ismth+1] = data[ii+ismth-buffer] * weights[ismth+1]
		end
		ndata[ii] = sum(smth)

	end

	return ndata

end

# ╔═╡ 4319fd0e-fd9f-424e-9286-3b3b5a844b73
function extract(geoname,iso,days)

	dystr  = @sprintf("%02d",days)

	ds = NCDataset(datadir(
		"wrf","processed",
		"$geoname-$(iso)QBUDGET-20190801_20201231.nc"
	))
	hvyp = smooth(dropdims(sum(reshape(ds["$(iso)P"][:],24,:),dims=1),dims=1),days)
	hvye = smooth(dropdims(sum(reshape(ds["$(iso)E"][:],24,:),dims=1),dims=1),days)
	close(ds)

	ds = NCDataset(datadir(
		"wrf","processed",
		"$geoname-QBUDGET-20190801_20201231.nc"
	))
	prcp = smooth(dropdims(sum(reshape(ds["P"][:],24,:),dims=1),dims=1),days)
	evap = smooth(dropdims(sum(reshape(ds["E"][:],24,:),dims=1),dims=1),days)
	close(ds)

	ds = NCDataset(datadir(
		"wrf","processed",
		"$geoname-$(iso)∇decompose-20190801_20201231.nc"
	))
	hvya = smooth(dropdims(mean(reshape(ds["$(iso)ADV"][:],24,:),dims=1),dims=1),days) * 86400
	close(ds)

	ds = NCDataset(datadir(
		"wrf","processed",
		"$geoname-∇decompose-20190801_20201231.nc"
	))
	advc = smooth(dropdims(mean(reshape(ds["ADV"][:],24,:),dims=1),dims=1),days) * 86400
	divg = smooth(dropdims(mean(reshape(ds["DIV"][:],24,:),dims=1),dims=1),days) * 86400
	close(ds)

	dsp = NCDataset(datadir(
		"wrf","processed",
		"$geoname-p_wwgt3-daily-20190801_20201231-smooth_$(dystr)days.nc"
	))
	pwgt = dsp["p_wwgt"][:] / 100
	pbl = dsp["P"][end,:] / 100 * 0.8; pbl[pbl.>800] .= 800
	close(dsp)

	return pbl .+ pwgt,prcp,evap,advc,hvyp,hvye,hvya,divg
	
end

# ╔═╡ 5bf90248-6ad6-4851-9c56-613d69f83d4b
function binning!(
	binstn,numstn,prcstn,
	rpnt, ppnt;
	ID, days=0, iso = "HDO", prfx = ""
)

	iso = "$(iso)_"
	stnstr = @sprintf("%02d",ID)
	
	geoname = "OTREC_wrf_$(prfx)$stnstr"
	pwgt,prcp,evap,advc,hvyp,hvye,hvya,divg = extract(geoname,iso,days)
	nt = length(prcp)

	for it = 1 : nt
		if (prcp[it]+advc[it]-evap[it]>2.5) && !isnan(pwgt[it])
			rind = argmin(abs.(prcp[it]+advc[it]-evap[it].-rpnt))
			pind = argmin(abs.(pwgt[it].-ppnt))
			numstn[rind,pind] += 1
			binstn[rind,pind] += hvyp[it] + hvya[it] - hvye[it]
			prcstn[rind,pind] += prcp[it] + advc[it] - evap[it]
		end
	end

	return

end

# ╔═╡ 9d38e14e-7226-4d57-ba6f-3b3382dfce1c
function threshold!(
	bin,num,prc;
	numthresh = 10
)

	bin[num.<numthresh] .= 0
	prc[num.<numthresh] .= 0

	return
	
end

# ╔═╡ 6fc8d69a-81d1-47c4-8609-8ec7914bc935
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
	ix = axesnum[2*ii].panel("l",width="0.6em",space=0)
	ix.pcolormesh(
		[0,1],wgtpbin,(tmpbin./tmpprc.-1)' * 1000,
		cmap="viridis",levels=lvls,extend="both"
	)
	ix.format(xlocator=[])

	ix = axesnum[2*ii-1].panel("r",width="0.6em",space=0)
	ix.pcolormesh(
		[0,1],wgtpbin,tmpnum',
		cmap="fire",levels=0:5:50,extend="both"
	)
	ix.format(xlocator=[])

	tmpbin = sum(iibin,dims=2)
	tmpprc = sum(iiprc,dims=2)
	tmpnum = sum(iinum,dims=2)
	threshold!(tmpbin,tmpnum,tmpprc)
	ix = axesnum[2*ii].panel("t",width="0.6em",space=0)
	ix.pcolormesh(
		rainbin,[0,1],(tmpbin./tmpprc.-1)' * 1000,
		cmap="viridis",levels=lvls,extend="both"
	)
	ix.format(ylocator=[])

	ix = axesnum[2*ii-1].panel("t",width="0.6em",space=0)
	ix.pcolormesh(
		rainbin,[0,1],tmpnum',
		cmap="fire",levels=0:5:50,extend="both"
	)
	ix.format(ylocator=[])

	threshold!(iibin,iinum,iiprc)
	c1 = axesnum[2*ii].pcolormesh(
		rainbin,wgtpbin,(iibin./iiprc .- 1)' * 1000,
		cmap="viridis",levels=lvls,extend="both"
	)
	c2 = axesnum[2*ii-1].pcolormesh(
		rainbin,wgtpbin,iinum',
		cmap="fire",levels=0:5:50,extend="both"
	)

	if returncinfo
		return c1,c2
	else
		return
	end

end

# ╔═╡ 987ab0fb-a376-4470-bcad-0e5681c6ca84
function axesformat!(axesnum)

	for ax in axesnum
		ax.format(
			ylim=(1000,0),xlim=(0,100),xlocator=0:100:200,ylocator=0:250:1000,
			xlabel=L"$P - E + u\cdot\nabla q$ / kg m$^{-2}$ day$^{-1}$",
			leftlabels=["ITCZ","CrossITCZ","PAC2ATL"]
		)
	end

	for ii = 1 : 5
		axesnum[2*ii].format(ultitle="$((ii-1)*5+1)-$(5*ii)")
		axesnum[2*ii+10].format(ultitle="$((ii-1)*5+1)-$(5*ii)")
		axesnum[2*ii+20].format(ultitle="$((ii-1)*5+1)-$(5*ii)")
	end
	
	naxs = length(axesnum)
	# axesnum[2].format(ylabel=L"$p_{sfc} - p'_{q\omega}$ / hPa")
	axesnum[12].format(ylabel=L"$p_{q\omega,bl}$ / hPa")
	# axesnum[22].format(ylabel=L"$p_{sfc} - p'_{q\omega}$ / hPa")

	return

end

# ╔═╡ 2fd946e2-bf3e-406f-9a19-5aa72b5d1640
begin
	rbin =  0 : 10 : 250; rpnt = (rbin[1:(end-1)] .+ rbin[2:end]) / 2
	pbin = 0 : 50 : 1000; ppnt = (pbin[1:(end-1)] .+ pbin[2:end]) / 2
	nr = length(rpnt); np = length(ppnt)
	abin = zeros(nr,np); anum = zeros(nr,np); aprc = zeros(nr,np)
	bbin = zeros(nr,np); bnum = zeros(nr,np); bprc = zeros(nr,np)
	cbin = zeros(nr,np); cnum = zeros(nr,np); cprc = zeros(nr,np)
	dbin = zeros(nr,np); dnum = zeros(nr,np); dprc = zeros(nr,np)
	ebin = zeros(nr,np); enum = zeros(nr,np); eprc = zeros(nr,np)
	fbin = zeros(nr,np); fnum = zeros(nr,np); fprc = zeros(nr,np)
	gbin = zeros(nr,np); gnum = zeros(nr,np); gprc = zeros(nr,np)
	hbin = zeros(nr,np); hnum = zeros(nr,np); hprc = zeros(nr,np)
	ibin = zeros(nr,np); inum = zeros(nr,np); iprc = zeros(nr,np)
	jbin = zeros(nr,np); jnum = zeros(nr,np); jprc = zeros(nr,np)
	kbin = zeros(nr,np); knum = zeros(nr,np); kprc = zeros(nr,np)
	lbin = zeros(nr,np); lnum = zeros(nr,np); lprc = zeros(nr,np)
	mbin = zeros(nr,np); mnum = zeros(nr,np); mprc = zeros(nr,np)
	nbin = zeros(nr,np); nnum = zeros(nr,np); nprc = zeros(nr,np)
	obin = zeros(nr,np); onum = zeros(nr,np); oprc = zeros(nr,np)
	md"Preallocation of arrays ..."
end

# ╔═╡ c57ae725-3056-481c-a004-a916192744be
# ╠═╡ show_logs = false
begin
	abin .= 0; aprc .= 0; anum .= 0
	bbin .= 0; bprc .= 0; bnum .= 0
	cbin .= 0; cprc .= 0; cnum .= 0
	dbin .= 0; dprc .= 0; dnum .= 0
	ebin .= 0; eprc .= 0; enum .= 0
	fbin .= 0; fprc .= 0; fnum .= 0
	gbin .= 0; gprc .= 0; gnum .= 0
	hbin .= 0; hprc .= 0; hnum .= 0
	ibin .= 0; iprc .= 0; inum .= 0
	jbin .= 0; jprc .= 0; jnum .= 0
	kbin .= 0; kprc .= 0; knum .= 0
	lbin .= 0; lprc .= 0; lnum .= 0
	mbin .= 0; mprc .= 0; mnum .= 0
	nbin .= 0; nprc .= 0; nnum .= 0
	obin .= 0; oprc .= 0; onum .= 0
	
	for istn = 1 : 5
		binning!(abin,anum,aprc,rpnt,ppnt,ID=istn,prfx="ITCZ",iso="O18",days=7)
	end
	for istn = 6 : 10
		binning!(bbin,bnum,bprc,rpnt,ppnt,ID=istn,prfx="ITCZ",iso="O18",days=7)
	end
	for istn = 11 : 15
		binning!(cbin,cnum,cprc,rpnt,ppnt,ID=istn,prfx="ITCZ",iso="O18",days=7)
	end
	for istn = 16 : 20
		binning!(dbin,dnum,dprc,rpnt,ppnt,ID=istn,prfx="ITCZ",iso="O18",days=7)
	end
	for istn = 21 : 25
		binning!(ebin,enum,eprc,rpnt,ppnt,ID=istn,prfx="ITCZ",iso="O18",days=7)
	end
	for istn = 1 : 5
		binning!(fbin,fnum,fprc,rpnt,ppnt,ID=istn,prfx="CrossITCZ",iso="O18",days=7)
	end
	for istn = 6 : 10
		binning!(gbin,gnum,gprc,rpnt,ppnt,ID=istn,prfx="CrossITCZ",iso="O18",days=7)
	end
	for istn = 11 : 15
		binning!(hbin,hnum,hprc,rpnt,ppnt,ID=istn,prfx="CrossITCZ",iso="O18",days=7)
	end
	for istn = 16 : 20
		binning!(ibin,inum,iprc,rpnt,ppnt,ID=istn,prfx="CrossITCZ",iso="O18",days=7)
	end
	for istn = 21 : 25
		binning!(jbin,jnum,jprc,rpnt,ppnt,ID=istn,prfx="CrossITCZ",iso="O18",days=7)
	end
	for istn = 1 : 5
		binning!(kbin,knum,kprc,rpnt,ppnt,ID=istn,prfx="PAC2ATL",iso="O18",days=7)
	end
	for istn = 6 : 10
		binning!(lbin,lnum,lprc,rpnt,ppnt,ID=istn,prfx="PAC2ATL",iso="O18",days=7)
	end
	for istn = 11 : 14
		binning!(mbin,mnum,mprc,rpnt,ppnt,ID=istn,prfx="PAC2ATL",iso="O18",days=7)
	end
	for istn = 16 : 20
		binning!(nbin,nnum,nprc,rpnt,ppnt,ID=istn,prfx="PAC2ATL",iso="O18",days=7)
	end
	for istn = 21 : 25
		binning!(obin,onum,oprc,rpnt,ppnt,ID=istn,prfx="PAC2ATL",iso="O18",days=7)
	end
	
	pplt.close(); f2,a2 = pplt.subplots(
		[[2,1,4,3,6,5,8,7,10,9],
		[12,11,14,13,16,15,18,17,20,19],
		[22,21,24,23,26,25,28,27,30,29]],
		aspect=0.5,axwidth=0.5,wspace=[0,1.5,0,1.5,0,1.5,0,1.5,0]
	)

	c2_1,c2_2 = 
	plotbin!(a2,01,rbin,pbin,abin,aprc,anum,-20:-8,returncinfo=true)
	plotbin!(a2,02,rbin,pbin,bbin,bprc,bnum,-20:-8)
	plotbin!(a2,03,rbin,pbin,cbin,cprc,cnum,-20:-8)
	plotbin!(a2,04,rbin,pbin,dbin,dprc,dnum,-20:-8)
	plotbin!(a2,05,rbin,pbin,ebin,eprc,enum,-20:-8)
	plotbin!(a2,06,rbin,pbin,fbin,fprc,fnum,-20:-8)
	plotbin!(a2,07,rbin,pbin,gbin,gprc,gnum,-20:-8)
	plotbin!(a2,08,rbin,pbin,hbin,hprc,hnum,-20:-8)
	plotbin!(a2,09,rbin,pbin,ibin,iprc,inum,-20:-8)
	plotbin!(a2,10,rbin,pbin,jbin,jprc,jnum,-20:-8)
	plotbin!(a2,11,rbin,pbin,kbin,kprc,knum,-20:-8)
	plotbin!(a2,12,rbin,pbin,lbin,lprc,lnum,-20:-8)
	plotbin!(a2,13,rbin,pbin,mbin,mprc,mnum,-20:-8)
	plotbin!(a2,14,rbin,pbin,nbin,nprc,nnum,-20:-8)
	plotbin!(a2,15,rbin,pbin,obin,oprc,onum,-20:-8)

	axesformat!(a2)
	a2[1].format(suptitle="7-Day WRF Moving Average")

	f2.colorbar(c2_1,row=[1],locator=-20:4:-8,label=L"$\delta^{18}$O / $\perthousand$",minorlocator=-150:5:-45)
	f2.colorbar(c2_2,row=[2],locator=0:10:50,label="Number of Observations")
	
	f2.savefig(projectdir("figures","fig8-pqblomega.png"),transparent=false,dpi=400)
	load(projectdir("figures","fig8-pqblomega.png"))
end

# ╔═╡ Cell order:
# ╠═2e7c33da-f8b5-11ec-08f2-2581af96575f
# ╟─e32a00ee-5f32-47a1-a983-91fb77bc5d18
# ╟─bec4e6f2-c2ea-421e-8c57-33e1ef90aa21
# ╟─59c930cd-5b7f-4047-8660-615148d1bd9f
# ╟─441f47a7-5757-4b24-8b52-a2877e0f0287
# ╟─f1720645-69a8-4f45-a6c1-8c06279d3590
# ╟─4319fd0e-fd9f-424e-9286-3b3b5a844b73
# ╟─5bf90248-6ad6-4851-9c56-613d69f83d4b
# ╟─9d38e14e-7226-4d57-ba6f-3b3382dfce1c
# ╟─6fc8d69a-81d1-47c4-8609-8ec7914bc935
# ╟─987ab0fb-a376-4470-bcad-0e5681c6ca84
# ╟─2fd946e2-bf3e-406f-9a19-5aa72b5d1640
# ╟─c57ae725-3056-481c-a004-a916192744be

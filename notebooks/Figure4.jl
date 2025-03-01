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
	using GeoRegions
	using NCDatasets
	using Printf
	using Statistics
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("ultraplot")

	include(srcdir("common.jl"))
	
	md"Loading modules for the 2024GLConvIso project..."
end

# ╔═╡ 2e7c33da-f8b5-11ec-08f2-2581af96575f
md"
# 03b. Station W-Weighted Pressure, $\sigma$
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
		"$geoname-p_wwgt-daily-20190801_20201231-smooth_$(dystr)days.nc"
	))
	pwgt = dsp["p_wwgt"][:] / 100
	pwgt[(pwgt.>1000).|(pwgt.<0)] .= NaN
	pwgt = dsp["σ_wwgt"][:]
	pwgt[(pwgt.>1).|(pwgt.<0)] .= NaN
	close(dsp)

	
	# dsp = NCDataset(datadir(
	# 	"wrf3","processed",
	# 	"$geoname-wcoeff-daily-20190801_20201231-smooth_$(dystr)days.nc"
	# ))
	# wc = dsp["wcoeff"][:,:]
	# close(dsp)
	# wcconv = sum(wc[:,[1,3,5,7,9]],dims=2)[:];
	# wctpbt = sum(wc[:,[2,4,6,8,10]],dims=2)[:];
	# wcabs  = sqrt.(wcconv.^2 .+ wctpbt.^2)[:];
	# wctpbt = wctpbt ./ wcabs
	# ii = wcconv .< 0

	return pwgt,prcp,evap,advc,hvyp,hvye,hvya,divg
	
end

# ╔═╡ 5bf90248-6ad6-4851-9c56-613d69f83d4b
function binning!(
	binstn,numstn,prcstn,
	rpnt, ppnt;
	ID, days=0, iso = "HDO", box = false, bnum = 1
)

	iso = "$(iso)_"
	stnstr = @sprintf("%02d",ID)

	valid = readdlm(datadir("wrf","wrfvscalc.txt"))[ID,:]

	if box
		for ibox = 1 : bnum
			if valid[ibox] < 0.15
				boxstr = @sprintf("%0d",ibox)
				geoname = "OTREC_wrf_stn$(stnstr)_box$(boxstr)"
				pwgt,prcp,evap,advc,hvyp,hvye,hvya,divg = extract(geoname,iso,days)
				nt = length(prcp)
			
				for it = 1 : nt
					if ((prcp[it]+advc[it]-evap[it])>2.5) && !isnan(pwgt[it])# && (advc[it]>0) && (divg[it]<0)
						rind = argmin(abs.(prcp[it]+advc[it]-evap[it].-rpnt))
						pind = argmin(abs.(pwgt[it].-ppnt))
						numstn[rind,pind] += 1
						binstn[rind,pind] += hvyp[it] + hvya[it] - hvye[it]
						prcstn[rind,pind] += prcp[it] + advc[it] - evap[it]
					end
				end
			end
		end
	else
		geoname = "OTREC_wrf_stn$stnstr"
		pwgt,prcp,evap,advc,hvyp,hvye,hvya,divg = extract(geoname,iso,days)
		nt = length(prcp)
	
		for it = 1 : nt
			if (prcp[it]>2.5) && (advc[it]>0) && (divg[it]<0) && !isnan(pwgt[it])
				rind = argmin(abs.(prcp[it]+advc[it]-evap[it].-rpnt))
				pind = argmin(abs.(pwgt[it].-ppnt))
				numstn[rind,pind] += 1
				binstn[rind,pind] += hvyp[it] + hvya[it] - hvye[it]
				prcstn[rind,pind] += prcp[it] + advc[it] - evap[it]
			end
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

# ╔═╡ 1343fbae-0ebd-4237-8273-0ebab8325424
function axesformat!(axesnum)

	for ax in axesnum
		ax.format(
			ylim=(1,0),xlim=(0,100),xlocator=0:100:200,ylabel=L"$\sigma_\omega$",
			xlabel=L"$P - E + u\cdot\nabla q$ / kg m$^{-2}$ day$^{-1}$"
		)
	end

	axesnum[2].format(ultitle="(a) All Stations")
	axesnum[2].format(ultitle="(a) Colombia")
	axesnum[4].format(ultitle="(b) San Andres")
	axesnum[6].format(ultitle="(c) Buenaventura")
	axesnum[6].text(45,0.270,"Bahia Solano",fontsize=10)
	axesnum[8].format(ultitle="(d) Quibdo")
	axesnum[10].format(ultitle="(e) Costa Rica")
	axesnum[12].format(ultitle="(f) EEFMB")
	axesnum[12].text(39,0.270,"ADMQ",fontsize=10)
	axesnum[12].text(39,0.390,"CGFI",fontsize=10)
	axesnum[14].format(ultitle="(g) Cahuita")
	axesnum[14].text(45,0.270,"Bataan",fontsize=10)
	axesnum[14].text(45,0.390,"Limon",fontsize=10)
	axesnum[16].format(ultitle="(h) Liberia")
	axesnum[16].text(45,0.270,"OSA",fontsize=10)

	return

end

# ╔═╡ 2fd946e2-bf3e-406f-9a19-5aa72b5d1640
begin
	rbin = 0 : 10 : 250; rpnt = (rbin[1:(end-1)] .+ rbin[2:end]) / 2
	pbin = 0 : 0.050 : 1; ppnt = (pbin[1:(end-1)] .+ pbin[2:end]) / 2
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
	md"Preallocation of arrays ..."
end

# ╔═╡ b6500812-fd5e-4842-8855-655822d170f4
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
	
	for istn = 1
		binning!(bbin,bnum,bprc,rpnt,ppnt,ID=istn,iso="O18",box=true,bnum=4,days=7)
	end
	for istn = [3,4]
		binning!(cbin,cnum,cprc,rpnt,ppnt,ID=istn,iso="O18",box=true,bnum=4,days=7)
	end
	for istn = 2
		binning!(dbin,dnum,dprc,rpnt,ppnt,ID=istn,iso="O18",box=true,bnum=4,days=7)
	end
	for istn = 5 : 7
		binning!(ebin,enum,eprc,rpnt,ppnt,ID=istn,iso="O18",box=true,bnum=4,days=7)
	end
	for istn = 9 : 11
		binning!(fbin,fnum,fprc,rpnt,ppnt,ID=istn,iso="O18",box=true,bnum=4,days=7)
	end
	for istn = [8,12]
		binning!(gbin,gnum,gprc,rpnt,ppnt,ID=istn,iso="O18",box=true,bnum=4,days=7)
	end

	hbin .= bbin.+cbin.+dbin; ibin .= ebin.+fbin.+gbin
	hprc .= bprc.+cprc.+dprc; iprc .= eprc.+fprc.+gprc
	hnum .= bnum.+cnum.+dnum; inum .= enum.+fnum.+gnum
	
	pplt.close(); f1,a1 = pplt.subplots(
		[[2,1,4,3,6,5,8,7],[10,9,12,11,14,13,16,15]],
		aspect=0.5,axwidth=0.6,wspace=[0,1.5,0,1.5,0,1.5,0]
	)

	c1_1,c1_2 = 
	plotbin!(a1,1,rbin,pbin,hbin,hprc,hnum,-22:-10,returncinfo=true)
	plotbin!(a1,2,rbin,pbin,bbin,bprc,bnum,-22:-10)
	plotbin!(a1,3,rbin,pbin,cbin,cprc,cnum,-22:-10)
	plotbin!(a1,4,rbin,pbin,dbin,dprc,dnum,-22:-10)
	plotbin!(a1,5,rbin,pbin,ibin,iprc,inum,-22:-10)
	plotbin!(a1,6,rbin,pbin,ebin,eprc,enum,-22:-10)
	plotbin!(a1,7,rbin,pbin,fbin,fprc,fnum,-22:-10)
	plotbin!(a1,8,rbin,pbin,gbin,gprc,gnum,-22:-10)

	axesformat!(a1)
	a1[1].format(suptitle="7-Day WRF Moving Average")

	f1.colorbar(c1_1,loc="r",rows=1,locator=-22:2:-8,label=L"$\delta^{18}$O / $\perthousand$",minorlocator=-150:5:-45)
	f1.colorbar(c1_2,loc="r",rows=2,locator=0:10:50,label="Number of Observations")
	
	f1.savefig(projectdir("figures","fig4-wrfadvec.png"),transparent=false,dpi=400)
	load(projectdir("figures","fig4-wrfadvec.png"))
end

# ╔═╡ Cell order:
# ╟─2e7c33da-f8b5-11ec-08f2-2581af96575f
# ╟─e32a00ee-5f32-47a1-a983-91fb77bc5d18
# ╟─bec4e6f2-c2ea-421e-8c57-33e1ef90aa21
# ╟─59c930cd-5b7f-4047-8660-615148d1bd9f
# ╟─441f47a7-5757-4b24-8b52-a2877e0f0287
# ╟─f1720645-69a8-4f45-a6c1-8c06279d3590
# ╟─4319fd0e-fd9f-424e-9286-3b3b5a844b73
# ╟─5bf90248-6ad6-4851-9c56-613d69f83d4b
# ╟─9d38e14e-7226-4d57-ba6f-3b3382dfce1c
# ╟─6fc8d69a-81d1-47c4-8609-8ec7914bc935
# ╟─1343fbae-0ebd-4237-8273-0ebab8325424
# ╟─2fd946e2-bf3e-406f-9a19-5aa72b5d1640
# ╟─b6500812-fd5e-4842-8855-655822d170f4

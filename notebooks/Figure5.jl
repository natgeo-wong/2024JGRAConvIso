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
	using NCDatasets
	using Printf
	using StatsBase
	
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

# ╔═╡ 4319fd0e-fd9f-424e-9286-3b3b5a844b73
function extract(geoname,days)

	dystr  = "daily-20190801_20201231-smooth_$(@sprintf("%02d",days))days"

	ds = NCDataset(datadir("wrf","processed","$geoname-dhodq-$(dystr).nc"))
	p = ds["P"][:,:] ./ 100
	dhdq = ds["dHDOdH2O"][:,:]
	dodq = ds["dO18dH2O"][:,:]
	close(ds)

	return p,dhdq,dodq
	
end

# ╔═╡ 887b52b3-8ea4-444d-8a6d-96fb6acb37f3
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

# ╔═╡ f224f902-cd89-4f4a-aedb-60ebbc79f4f5
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

# ╔═╡ 95de152a-39ce-448a-a931-ba393f86b629
function extractbudget(geoname,days)

	dystr  = "daily-20190801_20201231-smooth_$(@sprintf("%02d",days))days"

	ds = NCDataset(datadir(
		"wrf","processed",
		"$geoname-QBUDGET-20190801_20201231.nc"
	))
	prcp = smooth(dropdims(sum(reshape(ds["P"][:],24,:),dims=1),dims=1),days)
	evap = smooth(dropdims(sum(reshape(ds["E"][:],24,:),dims=1),dims=1),days)
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
		"$geoname-p_wwgt-$dystr.nc"
	))
	pwgt = dsp["p_wwgt"][:] / 100
	pwgt[(pwgt.>1000).|(pwgt.<0)] .= NaN
	pwgt = dsp["σ_wwgt"][:]
	pwgt[(pwgt.>1).|(pwgt.<0)] .= NaN
	close(dsp)

	return pwgt,prcp,evap,advc,divg
	
end

# ╔═╡ 5bf90248-6ad6-4851-9c56-613d69f83d4b
function plotdqdp(
	axes,ii;
	ID, days=0, cinfo = false
)

	dhdqbin = 0 : 0.02 : 1.2
	dodqbin = 0.9 : 0.002 : 1.02
	pbin    = vcat(0 : 50 : 550, 600 : 50 : 1000)
	binHDO = zeros(length(dhdqbin)-1,length(pbin)-1)
	binO18 = zeros(length(dodqbin)-1,length(pbin)-1)
	dhdqbinplt = (dhdqbin[1:(end-1)] .+ dhdqbin[2:end])/2
	dodqbinplt = (dodqbin[1:(end-1)] .+ dodqbin[2:end])/2
	pbinplt = (pbin[1:(end-1)] .+ pbin[2:end])/2

	for stn in ID
		wvc = readdlm(datadir("wrf","wrfvscalc-smooth30.txt"))[stn,:]
		qvl = readdlm(datadir("wrf","qbudgetdiff-smooth30.txt"))[stn,:]
		stnstr = @sprintf("%02d",stn)
		if (wvc[1] < 0.15) .& (qvl[1] < 0.05)
			geoname = "OTREC_wrf_stn$stnstr"
			p,dhdq,dodq = extract(geoname,days)
			pwgt,prcp,evap,advc,divg = extractbudget(geoname,days)
			it = ((prcp.+advc.-evap).>2.5) .& (pwgt.<1) .& (pwgt.>0)
			binHDO[:,:] += fit(
				Histogram,(dhdq[:,it][:],p[:,it][:]),(dhdqbin,pbin)
			).weights
			binO18[:,:] += fit(
				Histogram,(dodq[:,it][:],p[:,it][:]),(dodqbin,pbin)
			).weights
		end
		for ibox = 1 : 4
			if (wvc[ibox+1] < 0.15) .& (qvl[ibox+1] < 0.05)
				boxstr = @sprintf("%d",ibox)
				geoname = "OTREC_wrf_stn$(stnstr)_box$(boxstr)"
				p,dhdq,dodq = extract(geoname,days)
				pwgt,prcp,evap,advc,divg = extractbudget(geoname,days)
				it = ((prcp.+advc.-evap).>2.5) .& (pwgt.<1) .& (pwgt.>0)
				binHDO[:,:] += fit(
					Histogram,(dhdq[:,it][:],p[:,it][:]),(dhdqbin,pbin)
				).weights * 0.25
				binO18[:,:] += fit(
					Histogram,(dodq[:,it][:],p[:,it][:]),(dodqbin,pbin)
				).weights * 0.25
			end
		end
	end
	lvls = vcat(0.1,0.5,1,2:2:10,15:5:40)
	c1 = axes[2*ii].pcolormesh(dhdqbinplt,pbinplt,(binHDO./sum(binHDO,dims=1)*100)',extend="both",levels=lvls)
	c2 = axes[2*ii-1].pcolormesh(dodqbinplt,pbinplt,(binO18./sum(binO18,dims=1)*100)',extend="both",levels=lvls)

	if cinfo
		return c1,c2
	else
		return
	end

end

# ╔═╡ 1343fbae-0ebd-4237-8273-0ebab8325424
function axesformat!(axes)

	for ax in axes
		ax.format(
			xlim=(-0.2,1.2),ylim=(1000,100),ylabel="Pressure / hPa",
			xlabel=L"$\partial_pq_h/\partial_pq$ / VSMOW",
			suptitle="7-Day Moving Average"
		)
	end

	for ii = 1 : 6
		axes[2*ii].format(xlim=(0.5,1.1),lltitle="HDO")
		axes[2*ii-1].format(xlim=(0.92,1.01),lltitle="O18")
	end

	axes[2].format(ultitle="(b) San Andres")
	axes[4].format(ultitle="(c) Buenaventura")
	axes[4].text(0.715,300,"Bahia Solano",fontsize=10)
	axes[6].format(ultitle="(d) Quibdo")
	axes[8].format(ultitle="(f) EEFMB")
	axes[8].text(0.69,300,"ADMQ",fontsize=10)
	axes[8].text(0.69,390,"CGFI",fontsize=10)
	axes[10].format(ultitle="(g) Cahuita")
	axes[10].text(0.72,300,"Bataan",fontsize=10)
	axes[10].text(0.72,390,"Limon",fontsize=10)
	axes[12].format(ultitle="(h) Liberia")
	axes[12].text(0.72,300,"OSA",fontsize=10)

	return

end

# ╔═╡ 8c211620-d632-4f23-85f5-a702faf82270
# ╠═╡ show_logs = false
begin
	pplt.close(); fig,axs = pplt.subplots(
		[[2,1,4,3,6,5],[8,7,10,9,12,11]],aspect=0.5,axwidth=0.75,
		wspace=[0,2,0,2,0]
	)

	c1,_ =
	plotdqdp(axs,1,ID=1,days=7,cinfo=true)
	plotdqdp(axs,2,ID=3:4,days=7)
	plotdqdp(axs,3,ID=2,days=7)
	plotdqdp(axs,4,ID=5:7,days=7)
	plotdqdp(axs,5,ID=9:11,days=7)
	plotdqdp(axs,6,ID=[8,12],days=7)
	
	axesformat!(axs)

	fig.colorbar(c1,length=0.75,locator=vcat(0.1,0.5,1,2:2:10,15:5:40),label="Percentage of Observations per Pressure Bin / %")
	fig.savefig(projectdir("figures","fig5-dhdq.png"),transparent=false,dpi=400)
	load(projectdir("figures","fig5-dhdq.png"))
end

# ╔═╡ Cell order:
# ╟─2e7c33da-f8b5-11ec-08f2-2581af96575f
# ╟─e32a00ee-5f32-47a1-a983-91fb77bc5d18
# ╟─bec4e6f2-c2ea-421e-8c57-33e1ef90aa21
# ╟─59c930cd-5b7f-4047-8660-615148d1bd9f
# ╟─4319fd0e-fd9f-424e-9286-3b3b5a844b73
# ╟─95de152a-39ce-448a-a931-ba393f86b629
# ╟─887b52b3-8ea4-444d-8a6d-96fb6acb37f3
# ╟─f224f902-cd89-4f4a-aedb-60ebbc79f4f5
# ╠═5bf90248-6ad6-4851-9c56-613d69f83d4b
# ╟─1343fbae-0ebd-4237-8273-0ebab8325424
# ╠═8c211620-d632-4f23-85f5-a702faf82270

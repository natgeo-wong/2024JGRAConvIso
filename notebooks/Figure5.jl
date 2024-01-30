### A Pluto.jl notebook ###
# v0.19.37

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
	pplt = pyimport("proplot")

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
function extract(geoname,iso,days)

	dystr  = "-smooth_$(@sprintf("%02d",days))days"

	ds = NCDataset(datadir(
		"wrf","processed",
		"$geoname-$(iso)dhqdp$(dystr).nc"
	))
	p = ds["P"][:]
	hdq = ds["$(iso)hq"][:,:]
	dhqdp = ds["$(iso)dhqdp"][:,:]
	close(ds)

	return p,hdq,dhqdp * 1e2
	
end

# ╔═╡ 5bf90248-6ad6-4851-9c56-613d69f83d4b
function plotdqdp(
	axes,ii,dqdpbin;
	ID, days=0, box = false, bnum = 1, cinfo = false
)

	binHDO = zeros(length(dqdpbin)-1,50)
	binO18 = zeros(length(dqdpbin)-1,50)
	binplt = (dqdpbin[1:(end-1)] .+ dqdpbin[2:end])/2

	for stn in ID
		stnstr = @sprintf("%02d",stn)
		if box
			for ibox = 1 : bnum
				boxstr = @sprintf("%02d",ibox)
				geoname = "OTREC_STN$(stnstr)_$(boxstr)"
				p,hdq,dhqdp = extract(geoname,"HDO_",days)
				for jj = 1 : 50
					binHDO[:,jj] += fit(Histogram,dhqdp[jj,:],dqdpbin).weights
				end
				p,hdq,dhqdp = extract(geoname,"O18_",days)
				for jj = 1 : 50
					binO18[:,jj] += fit(Histogram,dhqdp[jj,:],dqdpbin./10).weights
				end
			end
		else
			geoname = "OTREC_STN$stnstr"
			p,hdq,dhqdp = extract(geoname,"HDO_",days)
			for jj = 1 : 50
				binHDO[:,jj] += fit(Histogram,dhqdp[jj,:],dqdpbin).weights
			end
			p,hdq,dhqdp = extract(geoname,"O18_",days)
			for jj = 1 : 50
				binO18[:,jj] += fit(Histogram,dhqdp[jj,:],dqdpbin./10).weights
			end
		end
	end
	lvls = [0,1,sqrt(2),2,2*sqrt(2.5),5,7*sqrt(2),10,10*sqrt(2),20,20*sqrt(2.5),50]
	c1 = axes[2*ii-1].pcolormesh(binplt,1:50,(binHDO./sum(binHDO,dims=1))'*100,extend="both",levels=lvls)
	c2 = axes[2*ii].pcolormesh(binplt./10,1:50,(binO18./sum(binO18,dims=1))'*100,extend="both",levels=lvls)

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
			ylim=(1,50),ylabel="Model Level",
			xlabel=L"$\partial_p(q_h/q)$ / $\perthousand$ hPa$^{-1}$",
			suptitle="7-Day Moving Average"
		)
	end

	for ii = 1 : 6
		axes[2*ii-1].format(xlim=(-0.2,1.2),lrtitle="HDO",xlocator=-0.5:0.5:1.5)
		axes[2*ii].format(xlim=(-0.2,1.2)./10,lrtitle="O18",xlocator=-0.05:0.05:.15)
	end

	axes[2].format(ultitle="(a) San Andres")
	axes[4].format(ultitle="(b) Buenaventura")
	axes[4].text(0.032,39.5,"Bahia Solano",fontsize=10)
	axes[6].format(ultitle="(c) Quibdo")
	axes[8].format(ultitle="(d) EEFMB")
	axes[8].text(0.032,39.5,"ADMQ",fontsize=10)
	axes[8].text(0.032,34.5,"CGFI",fontsize=10)
	axes[10].format(ultitle="(e) Cahuita")
	axes[10].text(0.032,39.5,"Bataan",fontsize=10)
	axes[10].text(0.032,34.5,"Limon",fontsize=10)
	axes[12].format(ultitle="(f) Liberia")
	axes[12].text(0.024,39.5,"OSA",fontsize=10)

	return

end

# ╔═╡ 8c211620-d632-4f23-85f5-a702faf82270
begin
	pplt.close(); fig,axs = pplt.subplots(
		[[2,1,4,3,6,5],[8,7,10,9,12,11]],aspect=0.5,axwidth=0.75,
		wspace=[0,1.5,0,1.5,0]
	)

	c1,_ =
	plotdqdp(axs,1,-1:0.05:1.5,ID=1,box=true,bnum=4,days=7,cinfo=true)
	plotdqdp(axs,2,-1:0.05:1.5,ID=3:4,box=true,bnum=4,days=7)
	plotdqdp(axs,3,-1:0.05:1.5,ID=2,box=true,bnum=4,days=7)
	plotdqdp(axs,4,-1:0.05:1.5,ID=5:7,box=true,bnum=4,days=7)
	plotdqdp(axs,5,-1:0.05:1.5,ID=9:11,box=true,bnum=4,days=7)
	plotdqdp(axs,6,-1:0.05:1.5,ID=[8,12],box=true,bnum=4,days=7)
	
	axesformat!(axs)

	fig.colorbar(c1,length=0.75,locator=[0,1,2,5,10,20,50],label="Probability / %")
	fig.savefig(projectdir("figures","fig5-dhqdp.png"),transparent=false,dpi=400)
	load(projectdir("figures","fig5-dhqdp.png"))
end

# ╔═╡ Cell order:
# ╟─2e7c33da-f8b5-11ec-08f2-2581af96575f
# ╟─e32a00ee-5f32-47a1-a983-91fb77bc5d18
# ╟─bec4e6f2-c2ea-421e-8c57-33e1ef90aa21
# ╟─59c930cd-5b7f-4047-8660-615148d1bd9f
# ╟─4319fd0e-fd9f-424e-9286-3b3b5a844b73
# ╟─5bf90248-6ad6-4851-9c56-613d69f83d4b
# ╟─1343fbae-0ebd-4237-8273-0ebab8325424
# ╟─8c211620-d632-4f23-85f5-a702faf82270

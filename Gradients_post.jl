using CSV
using DataFrames
using Statistics
using Gadfly
using Cairo
using StatsBase
using CairoMakie
using Makie, FileIO

#Set working directory
cd(raw"\\hydro-nas\Team\Projects\5630_DSC\Bacon Island Model\Model\20220121")

#Directory of gradients csv files
Grad_dir=raw"\\hydro-nas\Team\Projects\5630_DSC\GIS\raster\Numpy"

#Let's import gradients data DataFrames

grads0 = CSV.File(joinpath(Grad_dir, "DRNS2018.csv"))|> DataFrame
grads0_WLR = CSV.File(joinpath(Grad_dir, "WLR_DRNS2018.csv"))|> DataFrame

#Let's import toedrains indexes
toedrains= CSV.File(raw"\\hydro-nas\Team\Projects\5630_DSC\GIS\vector\GW Model\25ft\Drains\ToeDrains_Index_Flopy.csv")|> DataFrame

toegrads=semijoin(grads0,toedrains,on=[:k,:i,:j])
toegrads_WLR=semijoin(grads0_WLR,toedrains,on=[:k,:i,:j])

for year in 2019:1:2070
    grads_dum=CSV.File(joinpath(Grad_dir, "DRNS"*string(year)*".csv"))|> DataFrame
    grads_WLR_dum=CSV.File(joinpath(Grad_dir, "WLR_DRNS"*string(year)*".csv"))|> DataFrame
    toegrads_dum=semijoin(grads_dum,toedrains,on=[:k,:i,:j])
    toegrads_WLR_dum=semijoin(grads_WLR_dum,toedrains,on=[:k,:i,:j])
    append!(toegrads,toegrads_dum)
    append!(toegrads_WLR,toegrads_WLR_dum)
end


toegrads_all=toegrads
#Condition for filter
cond(Grad::Float64) = Grad > 1

toegrads=filter(:Grad=>cond,toegrads)

toegrads_WLR_GT1=filter(:Grad=>cond,toegrads_WLR)


#Group by year
gdf = groupby(toegrads, :Year)
gdf_WLR = groupby(toegrads_WLR_GT1, :Year)

#Count the number of gradients greater than 1
GT1=combine(gdf, nrow => :count)
GT1[!,:Fraction]=GT1[:, :count]/3732*100
GT1[!,:Scenario].="BAU"
GT1_WLR=copy(GT1)
GT1_WLR[GT1_WLR.Year.>2018,:Fraction].=0.0
GT1_WLR[!,:Scenario].="WLR"

append!(GT1,GT1_WLR)

myplot=plot(GT1,x=:Year, y=:Fraction,Geom.line, color=:Scenario,
Guide.ylabel("Exit gradients greater \nthan critical value (%)"),
style(major_label_font_size=16pt,minor_label_font_size=14pt,
line_width=2pt,key_title_font_size=16pt,key_label_font_size=14pt));
draw(PNG("CriticalGradients.png", 10inch, 5.625inch), myplot)

#2D Density plot
plot(toegrads_all, y="Grad",x="Year", Geom.density2d, Guide.colorkey("Density"))

plot(toegrads, y="Grad",x="Year", Geom.density2d, Guide.colorkey("Density"))


plot(toegrads_WLR, y="Grad",x="Year", Geom.density2d, Guide.colorkey("Density"))

plot(toegrads_WLR, y="top",x="Year", Geom.density2d, Guide.colorkey("Density"))

#2070
cond(Year::Int64) = Year == 2070
toegrads_2070=filter(:Year=>cond,toegrads_all)
toegrads_2070_BAU=copy(toegrads_2070)
toegrads_2070[!,:Scenario].="BAU"
toegrads_WLR_2070=filter(:Year=>cond,toegrads_WLR)
toegrads_WLR_2070[!,:Scenario].="WLR"
append!(toegrads_2070,toegrads_WLR_2070)

Density=plot(toegrads_2070, x="Grad", color="Scenario", xintercept=[1.0], Geom.density, 
Geom.vline(color="red",style=[:dash]),
Guide.ylabel("Density"),
Guide.xlabel("Gradient"),
style(major_label_font_size=16pt,minor_label_font_size=14pt,
line_width=2pt,key_title_font_size=16pt,key_label_font_size=14pt),
)

draw(PNG(raw"\\hydro-nas\Team\Projects\5630_DSC\Report\Shorter_Version\Densities2070.png", 
10inch, 5inch), Density)

f = Figure()
Axis(f[1, 1], xlabel = "Hydraulic Exit Gradients", ylabel = "Cumulative Frequency",
title = "Empirical Cumulative Distribution\nFunctions for 2070")
WLR=ecdfplot!(toegrads_WLR_2070[!,:Grad]; color=:orange)
BAU=ecdfplot!(toegrads_2070_BAU[!,:Grad])
Legend(f[1, 2],
    [WLR, BAU],
    ["WLR", "BAU"])


f

#Let's plot SLR timeseries
SLR=CSV.File(raw"\\hydro-nas\Team\Projects\5630_DSC\Bacon Island Model\Model\Final\SLR.csv")|> DataFrame

SLR=SLR[:,[:"Year",:"2_ft"]]

#Let's subset by years less than or equal to 2070
cond(Year::Int64) = Year <= 2070
SLR=filter(:Year=>cond,SLR)

cond(Year::Int64) = Year >= 2020
SLR=filter(:Year=>cond,SLR)

p=plot(SLR,x=:"Year",y=:"2_ft", Geom.line,
color=[colorant"plum2"],
Guide.ylabel("Tidal Stage (ft NAVD88)"),
style(major_label_font_size=16pt,minor_label_font_size=14pt,
line_width=2pt,key_title_font_size=16pt,key_label_font_size=14pt))

draw(PNG(raw"\\hydro-nas\Team\Projects\5630_DSC\Report\Shorter_Version\TidalStage.png", 5inch, 5inch), p)

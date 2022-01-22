using CSV
using DataFrames
using Statistics
using Gadfly
using Cairo
using StatsBase

#Set working directory
cd(raw"\\hydro-nas\Team\Projects\5630_DSC\Bacon Island Model\Model\20220121")

#Directory of gradients csv files
Grad_dir=raw"\\hydro-nas\Team\Projects\5630_DSC\GIS\raster\Numpy"

#Let's import gradients data DataFrames

grads0 = CSV.File(joinpath(Grad_dir, "DRNS2018.csv"))|> DataFrame

#Let's import toedrains indexes
toedrains= CSV.File(raw"\\hydro-nas\Team\Projects\5630_DSC\GIS\vector\GW Model\25ft\Drains\ToeDrains_Index_Flopy.csv")|> DataFrame

toegrads=semijoin(grads0,toedrains,on=[:k,:i,:j])

for year in 2019:1:2070
    grads_dum=CSV.File(joinpath(Grad_dir, "DRNS"*string(year)*".csv"))|> DataFrame
    toegrads_dum=semijoin(grads_dum,toedrains,on=[:k,:i,:j])
    append!(toegrads,toegrads_dum)
end

#Condition for filter
cond(Grad::Float64) = Grad > 1

toegrads=filter(:Grad=>cond,toegrads)

#Group by year
gdf = groupby(toegrads, :Year)

#Count the number of gradients greater than 1
GT1=combine(gdf, nrow => :count)
GT1[!,:Fraction]=GT1[:, :count]/3732
myplot=plot(GT1,x=:Year, y=:Fraction,Geom.line)

#2D Density plot
plot(toegrads, y="Grad",x="Year", Geom.density2d, Guide.colorkey("Density"))




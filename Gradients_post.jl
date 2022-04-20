using CSV
using DataFrames
#using Statistics
using Gadfly
#using Cairo
#using StatsBase
#using CairoMakie
#using Makie, FileIO
using Colors

#Set working directory
cd(raw"\\hydro-nas\Team\Projects\5630_DSC\Paper\2022_02\Figures")

#Directory of gradients csv files
Grad_dir=raw"\\hydro-nas\Team\Projects\5630_DSC\GIS\CSV"

year_0=2018
year_f=2070

Scenarios=["BAU","WLR"]

#Let's import toedrain gradients data DataFrames

Toedrains_all=CSV.read(joinpath(Grad_dir, "BAU_TOEDRNS2018.csv"),
DataFrame)

col_names=Symbol.(names(Toedrains_all))
col_types=eltype.(eachcol(Toedrains_all))
named_tuple = (; zip(col_names, type[] for type in col_types )...)

Toedrains_all = DataFrame(named_tuple)
Toedrains_all[!,"Scenario"]=String[]

###Let's loop to create a dataframe with all the gradients

for year in year_0:year_f
    for scenario in Scenarios
        #Let's import yearly dataframe
        Toedrains_dum=CSV.read(joinpath(Grad_dir, scenario*"_TOEDRNS"*string(year)*".csv"),
        DataFrame)
        Toedrains_dum[!,"Scenario"].=scenario
        append!(Toedrains_all,Toedrains_dum)
    end
end

CSV.write(joinpath(Grad_dir,"TOEDRNS_all.csv", Toedrains_all)
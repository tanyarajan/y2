#=
International Macro and Trade
Assignment 2
Tanya Rajan
-------------------------------------
Description: This file uses data associated with Head
and Mayer (2014) to run the fixed-effects log-linear version
of the gravity regression on all years from 1948-2006.
=#

###################
# Preamble
###################
for package in ["LaTeXStrings", "StatFiles", "FixedEffectModels",
  "StatsModels", "RegressionTables"]
  Pkg.add(package)
end
using LaTeXStrings,Random, Distributions, Plots, FixedEffectModels,
      Latexify, DataFrames, StatFiles, StatsModels, RegressionTables #Load packages

###################
# Estimation
###################

# opening cleaned .dta file
filepath = "$(homedir())/Documents/git/psets/trade/pset2" # change to run code
cd(filepath)
df = DataFrame(load("data/clean.dta"))

# creating smaller dataset for initial Julia run
dfnon = dropmissing(df)
dfnon = dfnon[(dfnon[:year] .>= 2004) ,:]

# initial run on small dataset
output1 = reg(dfnon, @formula(lnflow ~ lndist + contig + comlang_off +
    fe(expyear) + fe(impyear)), Vcov.robust())

# estimation and timing
@time output = reg(df, @formula(lnflow ~ lndist + contig + comlang_off +
    fe(expyear) + fe(impyear)), Vcov.robust())

# writing to latex
regtable(output, renderSettings = latexOutput("tables/table3_Julia.tex"),
statisticformat="%0.5f")

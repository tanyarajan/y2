#=
International Macro and Trade
Assignment 1
Tanya Rajan
-------------------------------------
This code solves the DFS 1977 model and includes
functions to:
1. output a plot of the DFS schedules
2. solve the DFS model
3. calculate welfare
4. calculate trade volume.
It also includes a grid search method to find
a b-schedule with the same trade volume but
different gains from trade.
Note: grid search can be turned off to make the
code run faster
=#

###################
# Preamble
###################

for package in ["CSV","LaTeXStrings","Interpolations","PyCall" ]
  Pkg.add(package)
end
using CSV,DelimitedFiles,LaTeXStrings,Interpolations,PyCall, Random,
      Distributions, Plots, Latexify, DataFrames #Load packages


# opening txt files
filepath = "$(homedir())/Documents/y2/trade/pset1" # change to run code
cd(filepath)
a = readdlm("DFS1977_example_a.txt")
b = readdlm("DFS1977_example_b.txt")[:,1]

#= change this to false to switch off grid search method
to save computational time (code will plug in the)
optimized values from the grid search for the baseline model
from the problem set instead of searching) =#
gridsearch = false

###################
# Plotting Function
###################
# Plotting the A and B functions
function plotter(a,b)

  # foreign unit cost / home unit cost
  Az = a[:,1]./a[:,2]

  # calculating nu(z)
  nu = cumsum(b, dims=1)
  Bz = nu./(1 .- nu)

  # censoring upper graph values for Bz (since they tend to infinity)
  nufrac = nu ./ sum(b)
  upper = findfirst(x -> x > .95, nufrac)[1]

  # plotting the two lines
  plotout = hcat(Az[1:upper], Bz[1:upper])
  plotmat = plot(plotout,
            label = ["A(z)" "B(z)"], lw=2,
            legend=:topleft, xlabel="z", ylabel="ω")
  return plotmat
end

# run plotter function and output graph
plotter(a, b)
png(join([filepath,"tr_abgraph.png"], "/"))


###################
# Solver function
###################
function DFS1977solver(a::Array{Float64,2},b::Array{Float64,1},L::Array{Float64,1},g::Float64)
  # extract N
  Na = size(a)[1]
  Nb = size(b)[1]

  # generating the A function
  A = a[:,1]./a[:,2]

  #### Error checking ####
  if size(a)[2] != 2 error("Vector a must be a Nx2 vector")
  elseif size(a)[1] <=2 error("Vector a must have N>2 values")
  elseif any(x -> x < 0, a) error("Elements of vector a must be non-negative")
  elseif issorted(A, rev=true) == false error("A(z) is not monotone decreasing")
  elseif Na != Nb error("Dimensions of a and b do not match")
  elseif any(x -> x <= 0, b) error("Vector b must be strictly positive")
  elseif abs(sum(b)- 1)>.001 error("Vector b must sum to 1")
  elseif g <= 0 || g > 1 error("Scalar g is not in (0,1]")
  end

  #### Plotting ####
  plotter(a,b)

  #### Solving the model ####
  # intializing values
  Azxg = A .* g
  Azdg = A ./ g
  zstar = zeros(Na)
  zbar = zeros(Na)
  λ = zeros(Na)
  λstar = zeros(Na)
  bint = cumsum(b)

  # looping over wage values of A*g and finding the good z
  # which gives the same wage value for A/g
  for i in 1:Na
    # starting guess
     wage = Azxg[i]
     zstar[i] = i
     zbar_v = argmin(abs.(Azdg .- wage))
     zbar[i] = zbar_v
     λ[i] = bint[zbar_v]
     λstar[i] = 1 - bint[i]
    end

    # solving for omega bar
    ω = ((1 .- λstar)./(1 .- λ)).*(L[1]./L[2])

    # minimizing distance between wage guess and ω
    idx = argmin(abs.(ω - Azxg))

 return zstar[idx], zbar[idx], ω[idx]
end


# run function for different g values
Ltest = [1.0, 1.0] # labor values
solve1 = DFS1977solver(a, b, Ltest, 1.0) # g = 1
solve9 = DFS1977solver(a, b, Ltest, 0.9) # g = 0.9

# output table to latex
solveout = DataFrame(Parameter = ["z^*", "z-bar", "omega"],
                            NoCost = round.(collect(solve1), digits=3),
                            TradeCost = round.(collect(solve9), digits=3))
solveout["Parameter"] = latexstring.(solveout["Parameter"])
write(joinpath(filepath, "tr_solve.tex"), latexify(solveout; env=:table))

###################
#  Welfare Calcs
###################
function DFS1977welfare(a::Array{Float64,2},b::Array{Float64,1},L::Array{Float64,1},g::Float64)

  # solving for trade outcomes
  trade = DFS1977solver(a,b,L,g)
  zstar = convert(Int64, trade[1])
  zbar = convert(Int64, trade[2])
  ω = trade[3]
  Az = a[:,1] ./ a[:,2]
  effwage = ω ./ Az

  # if no trade in equilibrium
  if zstar == 1.0 && zbar == size(b)[1]
    return 0, 0, 0, 0, 0, 0, 0, 0
  end

  # wage normalization
  wstar = (1 ./ ω)

  # autarky welfare
  home_a = -sum(b.*log.(a[:,2]))
  fore_a = log(wstar) - sum(b.*log.(wstar.*a[:,1]))

  # trade welfare
  home_t = - b[1:zbar]'*log.(a[1:zbar,2]) -
              b[zbar+1:end]'*log.((wstar.*a[zbar+1:end,2])./g)
  fore_t = log(wstar) - b[zstar:end]'*log.(wstar.*a[zstar:end,1]) -
              b[1:zstar]'*log.(a[1:zstar,2]./g)

  # gains from trade
  home_gft = home_t - home_a
  fore_gft = fore_t - fore_a

  return home_a, home_t, home_gft, fore_a, fore_t, fore_gft
end

# solving for welfare for various levels of g
Ltest = [1.0, 1.0] # labor vector
btest = copy
out1 = DFS1977welfare(a,b,Ltest, 1.0) # g = 1
out9 = DFS1977welfare(a,b,Ltest, 0.9) # g = .9

# output table to latex
welout = DataFrame(Welfare = ["Home Autarky", "Home Trade", "Home GFT",
                            "Foreign Autarky", "Foreign Trade", "Foreign GFT"],
                            g1 = round.(collect(out1), digits=3),
                            g9 = round.(collect(out9), digits=3))
welout["Welfare"] = latexstring.(welout["Welfare"])
write(joinpath(filepath, "tr_welcalcs.tex"), latexify(welout; env=:table))


###################
#  Trade Volume
###################
function DFS1977volume(a::Array{Float64,2},b::Array{Float64,1},L::Array{Float64,1},g::Float64)

  # solve the DFS model
  out = DFS1977solver(a,b,Ltest, g)
  zstar = convert(Int64, out[1])
  zbar = convert(Int64, out[2])
  ω = out[3]
  wstar = 1 ./ ω

  # get welfare calculations
  gft = (DFS1977welfare(a,b,Ltest, g)[3],DFS1977welfare(a,b,Ltest, g)[6])

  # volume calculation
  if out[1] == 1.0 && out[2] == size(b)[1] # if no trade in equilibrium
    volume = 0
  else
    volume = (1/g).*(sum(b[1:zstar].*wstar.*L[1]) + sum(b[zbar:end].*L[2]))
  end
  return gft[1], gft[2], volume
end

###### Same Volume, Different GFT ######

# baseline model
gfix = 0.9
startvol=DFS1977volume(a,b,Ltest,gfix)
startsolve=DFS1977solver(a,b,Ltest,gfix)
start=DFS1977volume(a,b,Ltest,gfix)[3] # original b

# changing A schedule - uniform home tech progress
anew = copy(a)
anew[:,2] = anew[:,2].*.5
DFS1977volume(anew,b,Ltest,gfix)


# changing b by multiplying by a constant, c, over different cutpoints, p
pass = [0.001:0.001:.2;] # values of c to test
outpass = fill(1.0, length(pass), length(b)-1) # output vector

# grid search
if gridsearch == true
  for c = 1:length(pass)
    for p = 1:length(b)-1
      bpass = vcat(b[1:p].*pass[c],b[p+1:end]) #calculating new b schedules
      bpass = bpass/sum(bpass) # normalization
      outpass[c,p] = DFS1977volume(anew,bpass,Ltest,gfix)[3]
    end
  end

  # minimize the difference in trade volume from the baseline model
  diff = outpass .- start
  min1 = argmin(abs.(diff))
  outpass[min1]
else min1=[76,143] # pre-optimized quantities to save computational time
end

# calculating the new volume and GFT
bfin = vcat(b[1:min1[2]].*pass[min1[1]],b[min1[2]+1:end])
bfin = bfin/sum(bfin) # normalization
newbsolve = DFS1977solver(anew, bfin, Ltest, gfix)
newbvol = DFS1977volume(anew,bfin,Ltest,gfix)

# output table to latex
newout = vcat(collect(newbsolve), collect(newbvol))
baseout = vcat(collect(startsolve), collect(startvol))
bsolve = DataFrame(Parameter = ["z^*", "z-bar", "omega", "GFT Home", "GFT Foreign", "Volume"],
                            Baseline= round.(collect(baseout), digits=3),
                            New = round.(collect(newout), digits=3))
bsolve["Parameter"] = latexstring.(bsolve["Parameter"])
write(joinpath(filepath, "tr_adjustb.tex"), latexify(bsolve; env=:table))

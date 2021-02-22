## Function: Discrete approximation of Dornbusch, Fischer and Samuelson AER, 1977
#First column of a represents foreign, second column represents home
# Set up Julia packages
import Pkg
run(mkdir -p DFS77project) ##Create a folder for this project (https://docs.julialang.org/en/v1.0.0/manual/running-external-programs/)
Pkg.activate(“./DFS77project”) ##Add packages to this project (https://docs.julialang.org/en/v1/stdlib/Pkg/index.html)
for package in [“CSV”,“LaTeXStrings”,“Interpolations”,“PyCall”,“Plots”,“Latexify”]
  Pkg.add(package)
end
using CSV,DelimitedFiles,LaTeXStrings,Interpolations,PyCall,Random,DataFrames,Distributions,Plots,Latexify  #Load packages
# Set up CD
filepath = “$(homedir())/Documents/PhD Economics/Second year/First quarter/International Macroeconomics and Trade/Assignments/Assignment 1"
cd(filepath)
# Open the txt
a = readdlm("DFS1977_example_a.txt")
b = readdlm("DFS1977_example_b.txt")
## First part: Plot Figure 1
function plotter(a,b)
# Parameters
# Relative unit labor requirement
A = a[:,1]./a[:,2]
# Total number of goods
N = length(a)
# Definition of nu: Fraction of income spent (anywhere) on those goods in which the home country has comparative advantage
nu= cumsum(b,dims=1)
# Labor
L = ones(2)
# Representation of the demand side
B = (nu ./ (1 .- nu)).*(L[1] ./ L[2])
# Make plot
figure = plot(hcat(A[1:140], B[1:140]))
#title = "Figure 1", label = ["A(z)" "B(z,L*/L)"], xlabel="z", ylabel="ω")
return figure
end
plotter(a,b)
png(join([filepath,“graph1.png”], “/”))
## Second part: Solver function
function DFS1977solver(a::Array{Float64,2},b::Array{Float64,1},L::Array{Float64,1},g::Float64)
# Plotting
plotter(a,b)
# Size of a and b
N = size(a)[1]
    ## 1. Verify that a has a dimension N-by-2, N>2, and is non-negative
    # 1.1. We need two countries
    if size(a)[2] < 2
        error(“Too few countries”)
    elseif size(a)[2] > 2
        error(“Too many countries”)
    end
    # 1.2. We need N>2
    if size(a)[1] < 2
        error(“N cannot be less than 2”)
    elseif size(a)[1] == 2
        error(“We need N to be bigger than 2")
    end
    # 1.3. We need a > 0
    if a < 0
        error(“We need a to be positive”)
    elseif a == 0
        error(“We need a to be more than 0")
    end
    ## 2. Verify that A=a[:,1]./a[:,2] is monotone decreasing (equation 1 in DFS)
    if all(diff(a(:,1)./a(:,2))<=0) != 1
        error(“A is not monotone decreasing”)
    end
    ## 3. Verify that b is a vector of dimension N (the same length as A), strictly positive, and that sum(b)==1
    # 3.1. We need b to be a vector of dimension N
    if length(b) != length(a)
        error(“We need b to be a vector of dimension N”)
    end
    # 3.2. We need b to be positive
    if b < 0
        error(“b is negative and has to be positive”)
    elseif b == 0
        error(“b is zero and has to be positive”)
    end
    # 3.3. We need sum(b) == 1
    if sum(b) != 1
        error(“We need b to sum up to one”)
    end
    ## 4. Verify that g is scalar (0,1] (as assumed in DFS III.B)
    if length(g) > 1
        error(“We need g to be a scalar”)
    elseif g <= 0 || g > 1
        error(“We need g to be between zero and 1, 1 included”)
    end
# Equilibrium
A_foreign = max(A.*g)
A_home = max(A./g)
zbar_star = zeros[N]
zbar = zeros[N]
lambda = zeros[N]
lambda_star = zeros[N]
# To graph this, we are going to map omega from A_home and then retrieve zbar and then A_foreign to get zbar_star
for j in 1:N
    w = A_foreign[j]
    zstar[j] = j
    zbar_abs_diff = argmin(abs.(A_home .- w))
    zbar[j] = zbar_abs_diff
    lambda[j] = cumsum(b)[zbar_abs_diff]
    lambda_star[j] = 1 .- cumsum(b)[j]
    dist = argmin(abs.(omega .- A_foreign))
end
omega_bar = (1 .- lambda_star)./(1 .- lambda).*(L[1]./L[2])
return zbar_star[dist], zbar[dist], omega_bar[dist]
end
a_n = a[1:50,:]
b_n = b[1:50]
rand_1 = rand(150)
rand_2 = rand(150)
random_1 = hcat(rand_1,rand_2)
random_2 = rand(150) .- 0.5
output = DFS1977solver(a, b, L, .9)
## Third part: Calculate Welfare
function DFS1977calculatewelfare(a::Array{Float64,2},b::Array{Float64,1},L::Array{Float64,1},g::Float64)
# Output
solver_result = DFS1977solver(a, b, L, g)
zbar_star = Int(solver_result[1])
zbar = Int(solver_result[2])
omega_bar = solver_result[3]
w_efficient = omega ./ A
wage = 1./omega
# Autarky
home_autarky = - sum(b[1:end].*log.(a[1:end,2]))
foreign_autarky = log(wage)- sum(b[1:end].*log.(wage.*a[1:end,1]))
home_trade = - b[1:zbar_star].*log.(a[1:zbar_star,2]) - b[zbar_star:end].*log.(wage.*a[1:zbar_star,2])
foreign_trade = log(wage) - b[zbar_star:end].*log.(wage.*a[zbar_star:end,1]) - b[1:zbar_star].*log.(a[1:zbar_star,2]./g)
# Calculate the gains from trade
home_gain = home_trade - home_autarky
foreign_gain = foreign_trade - foreign_autarky
# Return zeros if there is no trade
if zbar_star == 1 && zbar == size(b)[1]
    return 0, 0, 0, 0, 0, 0
end
# Return
return home_autarky, foreign_autarky, home_trade, foreign_trade, home_gain, foreign_gain
end
# Labor
L = [1,1]
# Output
output_g1 = DFS1977calculatewelfare(a,b,L,1)
output_g09 = DFS1977calculatewelfare(a,b,L,0.9)
output_welfare = DataFrame(Var = [“Home Autarky”, “Foreign Autarky”, “Home Trade”, “Foreign Trade”, “Home Gain”, “Foreign Gain”], g1 = collect(output_g1),g09 = collect(output_g09))
write(joinpath(filepath, “welfare.tex”), “/”),latexify(welout; env=:table, latex = false)
## Fourth part: Trade volume

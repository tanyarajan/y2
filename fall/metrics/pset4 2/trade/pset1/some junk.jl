
lot(hcat(Azxg, Azdg))
zstarT = findfirst(x -> x <= 0, Atest.*gtest .- ωguessT)
zbarT = findfirst(x -> x <= 0, Atest./gtest .- ωguessT)




# testing solver fn
Ltest = [1.0, 1.0]
Atest = a[:,1]./a[:,2]
nutest = cumsum(b, dims=1)
Btest = nutest./(1 .- nutest)
nufractest = nutest ./ sum(b)
uppertest = findfirst(x -> x > .95, nufractest)[1]
pouttest = hcat(Atest[1:uppertest], Btest[1:uppertest])



gtest = .4

gtest = .9
ωoutT = maximum(Atest .* gtest)
ωdifT = 1
tolT = .1
counter = 1
while counter < 10000 && ωdifT > tolT
        ωguessT = ωoutT
        zstarT = findfirst(x -> x <= 0, Atest.*gtest .- ωguessT)
        #print(zstarT)
        zbarT = findfirst(x -> x <= 0, Atest./gtest .- ωguessT)
        λT = sum(b[1:zstarT])
        λstarT = sum(b[zstarT:end])
        ωoutT = ((1- λstarT)./(1-λT)).*(Ltest[1,]./Ltest[2,])
        #print(ωimpT)
        #ωoutT = findfirst(x -> x < 0, ωimpT .- ωguessT)
        ωdifT = abs(ωoutT - ωguessT)
        print(counter)
        counter += 1
end

findfirst(x -> x <= 0, Atest.*gtest .- ωoutT)





counter = 1
ωoutT = maximum(Atest ./ gtest)
ωguessT = ωoutT
zstarT = findfirst(x -> x <= 0, Atest.*gtest .- ωguessT)
#print(zstarT)
zbarT = findfirst(x -> x <= 0, Atest./gtest .- ωguessT)
λT = sum(b[1:zbarT])
λstarT = sum(b[zbarT:end])
ωoutT = ((1- λstarT)./(1-λT)).*(Ltest[1,]./Ltest[2,])
#print(ωimpT)
#ωoutT = findfirst(x -> x < 0, ωimpT .- ωguessT)
ωdifT = abs(ωoutT - ωguessT)
print(counter)
counter += 1

wageline = fill(ωoutT,150,1)
plot(hcat(Atest, Atest.*gtest, Atest./gtest, wageline))




#= testing plotter fn:
Atest = a[:,1]./a[:,2]
nutest = cumsum(b, dims=1)
Btest = nutest./(1 .- nutest)
nufractest = nutest ./ sum(b)
uppertest = findfirst(x -> x > .95, nufractest)[1]
pouttest = hcat(Atest[1:uppertest], Btest[1:uppertest])
=#

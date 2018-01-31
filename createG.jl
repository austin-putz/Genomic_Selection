#==============================================================================#
# createG.jl
#==============================================================================#

#==============================================================================#
# Read in genotypes
#==============================================================================#

# Author:    Austin Putz
# Created:   Sept 27, 2017
# Modified:  Sept 27, 2017
# License:   GPLv2

#==============================================================================#
# Read in genotypes
#==============================================================================#

# load in DataFrames package
using DataFrames

println("----------------------------------------")
println("Reading in G Matrix...")
println("----------------------------------------")

# read in genotypes
@time geno = readtable("current_geno_ASReml.mkr", separator=' ', header=false);

println("----------------------------------------")
println("Converting Genotypes to Array...")
println("----------------------------------------")

# convert to an array
@time M = Array(geno[:,2:end]);

# get size to check column means
println("Size of the M matrix\n", size(M), "\n")

println("----------------------------------------")
println("Calculating column means...")
println("----------------------------------------")

# check for anything put 0/1/2 genotypes
# need to write check

# get column means of M matrix
@time col_means = mean(M, 1);

println("----------------------------------------")
println("Calculating P, Z, and G...")
println("----------------------------------------")

# repeat column means (2*p) for P matrix
@time P = repmat(col_means, size(M)[1]);

# subtract off column means
@time Z = M - P;

# calculate G
@time G = Z*Z';

# get allele frequencies
@time p = col_means / 2

# get q
@time q = 1 - p

# set sum 2pq
@time sum2pq = 2*dot(p,q)

# get G
@time G = G / sum2pq

# calculate condition number of G
@time cond_G = cond(G)

# print condition
println("The condition of G is: ", cond_G, "\n")

# get size of G
println("Size of the G matrix", size(G), "\n")

println("----------------------------------------")
println("Checking diagonals of G...")
println("----------------------------------------")

# get diags
Gdiags = diag(G)

# load UnicodePlots
using UnicodePlots

# histogram of my diagonals
println(histogram(Gdiags, bins=10, title="G Diagonals"))

println("----------------------------------------")
println("Checking allele frequencies...")
println("----------------------------------------")

# convert to Array{Float64,1} instead of Array{Float64,2}
p_vec = reshape(p, size(p)[2])

# histogram of my allele frequencies
println(histogram(p_vec, bins=50, title="Allele Frequencies"))

# reshope G in to 1 column to plot
G1 = reshape(G, length(G))

# plot it
println(histogram(G1, bins=50, title="All G Elements"))

println("----------------------------------------")
println("Plot of G matrix...")
println("----------------------------------------")

# load PyPlots package: install with Pkg.add("PyPlots")
using PyPlot

# print matrix
println(PyPlot.matshow(G))

# save plot
savefig("Gmatrix.png")

println("----------------------------------------")
println("Write out G matrix...")
println("----------------------------------------")

# round G
G2 = round.(G, 4)

# write a space delimited file
writedlm("Gmatrix.txt", G2)

println("----------------------------------------")
println("Fill new rounded matrix for full and half sibs...")
println("----------------------------------------")

# get rounded G matrix for full and half sibs
@time Ground = zeros(size(G)[1], size(G)[2])

for i in 1:size(Ground)[1]
	for j in 1:size(Ground)[2]

		# set diagonals to 1
		if i==j
			new_val = 1.0
		else 

			if G[i,j] < 0.15
				new_val = 0
			elseif G[i,j] >= 0.15 && G[i,j] < 0.40
				new_val = 0.25
			else G[i,j] >= 0.40
				new_val = 0.5
			end

		end

		# set new value to rounded value of G, to get back A matrix
		Ground[i, j] = new_val

	end
end

# plot "A" Matrix
println(PyPlot.matshow(Ground))

# save plot
savefig("Amatrix.png")

# write out "A" matrix, not really
writedlm("Amatrix.grm", Ground)



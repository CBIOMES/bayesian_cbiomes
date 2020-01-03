using PyPlot, Statistics, Distributions, Random

#################################################################

x = collect(-4:0.1:4);
pdfx = pdf.(Normal(0,1),x);
fig, ax = PyPlot.subplots()
ax.scatter(x,pdfx)
ax.set_ylabel("pdfx",fontsize=16)
ax.set_xlabel("x",fontsize=16);

#################################################################

sampx = rand(Normal(0,1),500);
fig, ax = PyPlot.subplots()
ax.scatter(collect(1:1:500),sampx)
ax.set_ylabel("sampx",fontsize=16);

#################################################################

fig, ax = PyPlot.subplots()
ax.hist(sampx,bins=10)
ax.set_ylabel("freq",fontsize=16)
ax.set_xlabel("sampx",fontsize=16);

#################################################################

using StatsPlots

#################################################################

StatsPlots.plot(density(sampx),xlabel = "beta1",ylabel = "Density")


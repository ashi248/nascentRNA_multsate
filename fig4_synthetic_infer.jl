## SSA result
using DelaySSAToolkit
using Random, Distributions
using StatsBase
using Interpolations
using Plots
using JLD2
using BlackBoxOptim
include("nascentRNA_function.jl")
include("TXModel.jl")

function Log_Likelihood_nascent(data,p)
      id = data[:,1] .+1
      N = maximum(id)
      p0 = [p;N]
      T = 1;
      prob = FSP_steady_infer(T,p0)
      prob = abs.(prob)
      tot_loss = -sum(log.(prob[id]))
      return (tot_loss)
end

##
# Figure 4 (A)(B)(C)

tt = collect(range(0,20000,20001))
tf = 20000.0
l1 = 3;l2 = 3; u = 4; v = 1; k1= 20; k2 = 0


# generate sample data using stochastic simulation
α1 = l2; μ1 = 1/v; α2 = l1; μ2 = 1/u
ρ = k1; α3 = 1.0; μ3 = 1.0; α4 = 1.0; μ4 = 1.0
parms = [α1, μ1, α2, μ2, α3, μ3, α4, μ4, ρ]
jsol = TX_model(parms, tf)

mRNA = jsol[3,:]
nodes = (jsol.t,)
mRNA_itp = Interpolations.interpolate(nodes,mRNA, Gridded(Constant{Previous}()))
mRNA0 = mRNA_itp(tt)
mRNA_expr = mRNA0[10001:20001]


# infere parameters using resample method
numb = 400
params = zeros(numb ,5)
for k = 1:numb
    print(k)
    mRNA_sample = sample(mRNA_expr[1:10000],10000, replace = true)
    df1 = Int.(round.(mRNA_sample))
    cost_function(p) = Log_Likelihood_nascent(df1, p)
    bound1 = Tuple{Float64, Float64}[(1, 5),(1, 5),(log(0.1), log(10)),(log(0.1), log(10)),(log(1), log(100))]
    result = bboptimize(cost_function;SearchRange = bound1,MaxSteps = 2000)
    results1 =  result.archive_output.best_candidate
    l1 = Int(round(results1[1]))
    l2 = Int(round(results1[2]))
    u = exp(results1[3])
    v = exp(results1[4])
    k1 = exp(results1[5])
    params[k,:] = [l1,l2,u,v,k1]
end

save_object("data/infer_nascent.jld2",params)


##
# figure 4 B
params2 = load_object("data/infer_nascent.jld2")

fs = sort(countmap(params2[:,1]))
fig1 = bar(collect(fs.keys),fs.vals ./sum(fs.vals),color = :green3, xlims = [0.5,5.5],
bar_edges = true, bar_width = 0.5,ylabel = "Probability", label = "",
xlabel = "L1",labelfontsize=15,tickfontsize = 12,legendfontsize=12,
grid=0,framestyle=:box,legend = (0.1,0.8))
vline!([3], lw = 5, color = "red",label= "True value")

fs = sort(countmap(params2[:,2]))
fig2 = bar(collect(fs.keys),fs.vals ./sum(fs.vals),color = :green3, xlims = [0.5,5.5],
bar_edges = true, bar_width = 0.5,ylabel = "Probability",label = "",
xlabel = "L2",labelfontsize=15,tickfontsize = 12,legendfontsize=12,
grid=0,framestyle=:box)
vline!([3],lw = 5, color = "red",label= "")

plot(fig1,fig2,layout = (2, 1),size = (500,400),
grid = false, fontfamily="Helvetica")

Plots.savefig("figure1/nascent_estimated_1.png")


##
# figure 4C
# error bar
params2 = load_object("data/infer_nascent.jld2")
p3 = params2[:,3]; p4 = params2[:,4]; p5 = params2[:,5]

a = ["u","v","k1"]
med_p = [median(p3);median(p4);median(p5)]
min_p = [quantile(p3,0.25);quantile(p4,0.25);quantile(p5,0.25)]
max_p = [quantile(p3,0.975);quantile(p4,0.975);quantile(p5,0.975)]

err1 = med_p .-min_p
err2 = max_p .-med_p

scatter(a, med_p,yerror=(err1,err2),msw = 2, ms=7,color = :green3,
yscale = :log,xlim = [0,3],label = "Median",
linewidth=2,size = (500,400),labelfontsize=18,tickfontsize = 16,legendfontsize=14,
grid=0,framestyle=:box,ylabel = "Estimated value")

scatter!([0.5,1.5,2.5],[4,1,20],markersize=8,markercolor = :red,
label = "True value",legend = (0.1,0.9))

Plots.savefig("figure1/nascent_estimated_2.png")

##
# figure 4A

p1 = params2[:,1];p2 = params2[:,2];
p3 = params2[:,3]; p4 = params2[:,4]; p5 = params2[:,5]
med_p = [median(p1),median(p2),median(p3),median(p4),median(p5)]

param_cyclic = med_p
mRNA_prob = proportionmap(mRNA_expr)
bar(mRNA_prob,label = "Synthetic data",color = "skyblue")

l1 = convert(Int,med_p[1])
l2 = convert(Int,med_p[2])
k1 = med_p[3]
u = med_p[4]
v = med_p[5]

prob1 = FSP_steady(l1,l2,1,k1,u,v,0,30)
N0 = 30
bins = collect(0:N0-1)
Plots.plot!(bins, prob1,size = (500,400),xlims = [0,20], dpi=600, linewidth = 5,
grid=0,framestyle=:box,label = "cyclic model",labelfontsize=18,tickfontsize = 14,
legendfontsize=14, color = :red,
xlabel = "Nascent RNA",ylabel = "Probability")


#-------------------------------------
mRNA_sample = mRNA_expr[1:10000]
df1 = Int.(round.(mRNA_sample))
cost_function(p) = Log_Likelihood_nascent(df1, p)
bound1 = Tuple{Float64, Float64}[(1, 1),(1, 1),(log(0.1), log(10)),(log(0.1), log(10)),(log(1), log(100))]
result = bboptimize(cost_function;SearchRange = bound1,MaxSteps = 10000)
results1 =  result.archive_output.best_candidate

l1 = Int(round(results1[1]))
l2 = Int(round(results1[2]))
u = exp(results1[3])
v = exp(results1[4])
k1 = exp(results1[5])

param_2state= [l1,l2,u,v,k1]

prob1 = FSP_steady(l1,l2,1,u,v,k1,0,30)
N0 = 30
bins = collect(0:N0-1)
Plots.plot!(bins, prob1,size = (500,400),dpi=600, linewidth = 5,
grid=0,framestyle=:box,label = "2-state model",labelfontsize=18,tickfontsize = 14,
legendfontsize=14, color = :blue,ls = :dash)

Plots.savefig("figure1/nascent_simulation_infer1.png")

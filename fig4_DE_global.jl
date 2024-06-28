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
using DataFrames
using StatsPlots
using Measures


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
#figure 4 (D)(E)

M0 = 600
U = rand(LogUniform(0.1,10),M0)
fc = sample(1:2:5,M0)
V =  U ./fc
K1 = rand(LogUniform(20,60),M0)
L1 =  sample(1:5,M0)
L2 =  sample(1:5,M0)

params_sample= DataFrame(u = U,v = V,k1= K1,l1 = L1, l2 = L2, fc = fc)

##
params_infer = zeros(M0,5)

for i=1:M0
    print(i)
    tt = collect(range(0,20000,20001))
    tf = 20000.0
    l1 = params_sample.l1[i];l2 = params_sample.l2[i]
    u = params_sample.u[i]; v = params_sample.v[i];
    k1= params_sample.k1[i]; k2 = 0

    # parameter values
    α1 = l2; μ1 = 1/v; α2 = l1; μ2 = 1/u
    ρ = k1; α3 = 1.0; μ3 = 1.0; α4 = 1.0; μ4 = 1.0

    parms = [α1, μ1, α2, μ2, α3, μ3, α4, μ4, ρ]

    jsol = TX_model(parms, tf)

    mRNA = jsol[3,:]
    nodes = (jsol.t,)
    mRNA_itp = Interpolations.interpolate(nodes,mRNA, Gridded(Constant{Previous}()))
    mRNA0 = mRNA_itp(tt)
    mRNA_expr = mRNA0[10001:20001]


    mRNA_sample = mRNA_expr
    df1 = Int.(round.(mRNA_sample))
    cost_function(p) = Log_Likelihood_nascent(df1, p)
    bound1 = Tuple{Float64, Float64}[(1, 5),(1, 5),(log(0.1), log(50)),(log(0.1), log(50)),(log(1), log(100))]
    result = bboptimize(cost_function;SearchRange = bound1, MaxSteps = 1000)
    results1 =  result.archive_output.best_candidate
    m = Int(round(results1[1]))
    n = Int(round(results1[2]))
    a1 = exp(results1[3])
    a2 = exp(results1[4])
    a3 = exp(results1[5])
    params_infer[i,:] = [a1,a2,a3,m,n]
end

rename!(params_sample,:l1 => :k1,:m => :L1,:n => :L2,:fon => :fc)
save_object("data/infer_nascent_global.jld2",params_sample)
save_object("data/infer_nascent_global_infer.jld2",params_infer)

##
params_sample = load_object("data/infer_nascent_global.jld2")
params_infer = load_object("data/infer_nascent_global_infer.jld2")

err_u = abs.(params_infer[:,1] .- params_sample.u)./params_sample.u
err_v = abs.(params_infer[:,2] .- params_sample.v)./params_sample.v
err_k1 = abs.(params_infer[:,3] .- params_sample.k1)./params_sample.k1
err_L1 = abs.(params_infer[:,4] .- params_sample.L1)./params_sample.L1
err_L2 = abs.(params_infer[:,5] .- params_sample.L2)./params_sample.L2


df_error = DataFrame(L1 = err_L1,L2 = err_L2,k1 = err_k1,
u = err_u, v = err_v,fc = params_sample.fc)
df_error_stack = stack(df_error,1:5)


error_gdf = groupby(df_error,:fc)
error_med = combine(error_gdf, [:L1,:L2,:k1,:u,:v].=>median,renamecols=false)

df_error_stack = stack(error_med,2:6)
@df df_error_stack  StatsPlots.groupedbar(:variable, :value,
 group = :fc,linewidth=2,outliers = false,
 labelfontsize=24,tickfontsize = 18,legendfontsize=18,
 ylabel = "Median Relative Error",label = ["fon=1/2" "fon=1/4" "fon=1/6"],
 size = (800, 500),legend = :topright,framestyle = :box,ylims = [0,0.6],
 left_margin = 8mm)

 Plots.savefig("figure1/nascent_simulation_medianRE_global.png")

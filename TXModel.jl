using DelaySSAToolkit
using Distributions
using Interpolations
using StatsBase

function TX_model(parms::Vector{Float64}, tf::Float64)
    # 1 OFF; 2 ON, 3 nascentRNA, 4 matureRNA

    # markovian part
    # ON --> ON + nascentRNA  Markovian reaction with rate ρ

    # non-markovian part
    # OFF => ON                 non-Markovian reaction, Gamma distribution with shape α1 and mean μ1
    # ON => OFF                 non-Markovian reaction, Gamma distribution with shape α2 and mean μ2
    # nascentRNA => matureRNA   non-Markovian reaction, Gamma distribution with shape α3 and mean μ3
    # matureRNA => 0            non-Markovian reaction, Gamma distribution with shape α4 and mean μ4

    # Gamma distribution referring https://juliastats.org/Distributions.jl/v0.14/univariate.html#Distributions.Gamma

    α1 = parms[1]
    μ1 = parms[2]
    α2 = parms[3]
    μ2 = parms[4]
    α3 = parms[5]
    μ3 = parms[6]
    α4 = parms[7]
    μ4 = parms[8]
    ρ  = parms[9]

    rates = [ρ]
    reactant_stoich = [[2=>1]]
    net_stoich = [[3=>1]]
    mass_jump = DelaySSAToolkit.MassActionJump(rates, reactant_stoich, net_stoich; scale_rates=false)
    jumpset = DelaySSAToolkit.JumpSet((), (), nothing, mass_jump)


    u0 = [1,0,0,0]  # [OFF, ON, nascentRNA, matureRNA]
    de_chan0 = [[1e-8], [], [], []]
    tspan = (0, tf)
    dprob = DiscreteProblem(u0, tspan)

    # nascentRNA => matureRNA trigger
    delay_trigger_affect1! = function (integrator, rng)
        α = α3
        μ = μ3
        τ = 1.0
        append!(integrator.de_chan[3], τ)
    end


    #  OFF =>  ON complete; ON => OFF trigger
    delay_complete_affect1! = function (integrator, rng)
        integrator.u[1] -= 1 # OFF state minus 1
        integrator.u[2] += 1 # ON state plus 1
        α = α2
        μ = μ2
        τ = rand(Gamma(α, μ / α))
        append!(integrator.de_chan[2], τ) # add to the delay channel
    end

    #  ON =>  OFF complete; OFF => ON trigger
    delay_complete_affect2! = function (integrator, rng)
        integrator.u[2] -= 1 # OFF state minus 1
        integrator.u[1] += 1 # ON state plus 1
        α = α1
        μ = μ1
        τ = rand(Gamma(α, μ / α))
        append!(integrator.de_chan[1], τ) # add to the delay channel
    end

    #  nascentRNA => matureRNA complete; matureRNA => 0 trigger
    delay_complete_affect3! = function (integrator, rng)
        integrator.u[3] -= 1 # nascentRNA minus 1
        integrator.u[4] += 1 # matureRNA plus 1
        α = α4
        μ = μ4
        τ = rand(Gamma(α, μ / α))
        append!(integrator.de_chan[4], τ) # add to the delay channel
    end

    #  matureRNA => 0 complete
    delay_complete_affect4! = function (integrator, rng)
        integrator.u[4] -= 1 # matureRNA minus 1
    end



    delay_trigger = Dict(1 => delay_trigger_affect1!)
    delay_complete = Dict(1 => delay_complete_affect1!, 2 => delay_complete_affect2!, 3 => delay_complete_affect3!,
        4 => delay_complete_affect4!)
    delay_interrupt = Dict()
    delayjumpset = DelayJumpSet(delay_trigger, delay_complete, delay_interrupt)


    djprob = DelayJumpProblem(dprob, DelayRejection(), jumpset, delayjumpset,
        de_chan0, save_positions=(true, true))
    sol1 = solve(djprob, DelaySSAToolkit.SSAStepper())
    sol1
end

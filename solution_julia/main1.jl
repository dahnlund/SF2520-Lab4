# Computer Exercise 4, David Ahnlund, Emil Gestsson

using Plots
using LinearAlgebra
include("hyperbolic.jl")
include("plot1b.jl")
include("plot1c1.jl")
include("plot1c2.jl")

a = 2; D = 10; T = 4
 
N = 100
lambda = 0.4

tau = 2

function u_b1(t)
    return sin.(2*pi*t/tau)
end


function u_b2(t)
    return sign.(sin.(2*pi*t/tau))
end

#Exact solutions:

function u_exact1(x,t)
    return (-u_b1(x'/a .- t) .* (x'.-a*t.<0))'
end
function u_exact2(x,t)
    return (-u_b2(x'/a .- t) .* (x'.-a*t.<0))'
end

#### 1B
# For 'g_sin'
plot1b(u_b1, a, N, lambda, D, T, u_exact1)


# For 'g_sq'
plot1b(u_b2, a, N, lambda, D, T, u_exact2)
####


#### 1C
plot1c1(0.05:0.05:0.5, u_b2, 500, 3);

plot1c2(0.4, u_b2, 100:100:1500, 3);

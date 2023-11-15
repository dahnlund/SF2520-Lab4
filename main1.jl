# Computer Exercise 4, David Ahnlund, Emil Gestsson

using Plots
using LinearAlgebra
include("hyperbolic.jl")

a = 2; D = 10; T = 4
 
N = 100
l = 0.4

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

# For 'g_sin'

t, x, u_lf1 = hyperbolic(1, a, u_b1, N, l, D, T)

_, _, u_lw1 = hyperbolic(2, a, u_b1, N, l, D, T)

_, _, u_up1 = hyperbolic(3, a, u_b1, N, l, D, T)

plot1 = plot(x,u_lf1[:,end], label = "Lax-Friedrichs")
plot!(x,u_lw1[:,end], label = "Lax-Wendroff")
plot!(x,u_up1[:,end], label = "Upwind")
plot!(x, u_exact1(x,t[end]), label = "Exact")
xlabel!("x")
title!("Plot of u(x,"*string(round(t[end],digits=4))*"), sin")
display(plot1)

# For 'g_sq'

_, _, u_lf2 = hyperbolic(1, a, u_b2, N, l, D, T)

_, _, u_lw2 = hyperbolic(2, a, u_b2, N, l, D, T)

_, _, u_up2 = hyperbolic(3, a, u_b2, N, l, D, T)

plot2 = plot(x,u_lf2[:,end], label = "Lax-Friedrichs")
plot!(x,u_lw2[:,end], label = "Lax-Wendroff")
plot!(x,u_up2[:,end], label = "Upwind")
plot!(x, u_exact2(x,t[end]), label = "Exact")
xlabel!("x")
title!("Plot of u(x,"*string(round(t[end],digits=4))*"), sq")
display(plot2)

"""
%%
figure
plot1c1(0.05:0.05:0.5, u_b2, 500, 3);

figure
plot1c2(0.4, u_b2, 100:100:1500, 3);
"""
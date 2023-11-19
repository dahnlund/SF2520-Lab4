using Plots

function plot1b(bound, a, N, lambda, D, T, exact_func)


t, x, u_lf1 = hyperbolic(1, a, bound, N, lambda, D, T)
_, _, u_lw1 = hyperbolic(2, a, bound, N, lambda, D, T)
_, _, u_up1 = hyperbolic(3, a, bound, N, lambda, D, T)

plot1 = plot(x,u_lf1[:,end], label = "Lax-Friedrichs", size = (800,600))
plot!(x,u_lw1[:,end], label = "Lax-Wendroff")
plot!(x,u_up1[:,end], label = "Upwind")
plot!(x, exact_func(x,t[end]), label = "Exact")
xlabel!("x")
title!("Plot of u(x,"*string(round(t[end],digits=4))*"), "*string(bound))
display(plot1)

end
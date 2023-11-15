# Computer Exercise 4, David Ahnlund, Emil Gestsson

using Plots

function plot1c2(CFL, bound, N, scheme)

    a = 2; D = 10; T = 4;
    plot_ = plot()
    for i in eachindex(N)

        _, x, u = hyperbolic(scheme,a,bound,N[i], CFL, D, T)

        plot!(x, u[:,end])
    end
    title!("Plot of u(x,T), N in [" * string(N[1])*","*string(N[end])*"]")
    xlabel!("x")
    display(plot_)

end
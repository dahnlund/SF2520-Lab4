# Computer Exercise 4, David Ahnlund, Emil Gestsson
using Plots

function plot1c1(CFL, bound, N, scheme)

    a = 2; D = 10; T = 4;
    plot_ = plot(size = (800,600))
    for i in eachindex(CFL)

        _, x, u = hyperbolic(scheme,a,bound,N, CFL[i], D, T)

        plot!(x, u[:,end])
    end
    title!("Plot of u(x,T), CFL in [" * string(CFL[1])*","*string(CFL[end])*"]")
    xlabel!("x")
    display(plot_)

end
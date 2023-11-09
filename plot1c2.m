%% Computer Exercise 4, David Ahnlund, Emil Gestsson

function plot1c2(CFL, bound, N, scheme)

a = 2; D = 10; T = 4;

for i = 1:length(N)

    [~, x, u] = hyperbolic1D(scheme,a,bound, N(i), CFL, D, T);

    plot(x, u(:,end)); hold on
end
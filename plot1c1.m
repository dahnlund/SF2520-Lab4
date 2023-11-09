%% Computer Exercise 4, David Ahnlund, Emil Gestsson

function plot1c1(CFL, bound, N, scheme)

a = 2; D = 10; T = 4;

for i = 1:length(CFL)

    [~, x, u] = hyperbolic1D(scheme,a,bound,N, CFL(i), D, T);

    plot(x, u(:,end)); hold on
end
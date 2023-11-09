%% Computer Exercise 4, David Ahnlund, Emil Gestsson
clc, clear variables;

a = 2; D = 10; T = 4;
 
N = 200; 
l = 0.2;

tau = 2;
u_b1 = @(t) sin(2*pi*t/tau);

u_b2 = @(t) sign(sin(2*pi*t/tau));


% For 'g_sin'

[t, x, u_lf1] = hyperbolic1D(1, a, u_b1, N, l, D, T);

[~, ~, u_lw1] = hyperbolic1D(2, a, u_b1, N, l, D, T);

[~, ~, u_up1] = hyperbolic1D(3, a, u_b1, N, l, D, T);

plot(x,u_lf1(:,end)); hold on; plot(x,u_lw1(:,end));plot(x,u_up1(:,end))
xlabel("x")
legend("Lax-Friedrichs","Lax-Wendroff", "Upwind")

% For 'g_sq'

[~, ~, u_lf2] = hyperbolic1D(1, a, u_b2, N, l, D, T);

[~, ~, u_lw2] = hyperbolic1D(2, a, u_b2, N, l, D, T);

[~, ~, u_up2] = hyperbolic1D(3, a, u_b2, N, l, D, T);

figure
plot(x,u_lf2(:,end)); hold on; plot(x,u_lw2(:,end));plot(x,u_up2(:,end))
xlabel("x")
legend("Lax-Friedrichs","Lax-Wendroff", "Upwind")









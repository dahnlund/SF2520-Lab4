%% Computer Exercise 4, David Ahnlund, Emil Gestsson
clc, clear variables;

a = 2; D = 10; T = 4;
 
N = 100; 
l = 0.4;

tau = 2;
u_b1 = @(t) sin(2*pi*t/tau);

u_b2 = @(t) sign(sin(2*pi*t/tau));

%Exact solutions:

u_exact1 = @(x,t) -u_b1(x'/a - t) .* (x'-a*t<0);
u_exact2 = @(x,t) -u_b2(x'/a - t) .* (x'-a*t<0);

% For 'g_sin'

[t, x, u_lf1] = hyperbolic1D(1, a, u_b1, N, l, D, T);

[~, ~, u_lw1] = hyperbolic1D(2, a, u_b1, N, l, D, T);

[~, ~, u_up1] = hyperbolic1D(3, a, u_b1, N, l, D, T);

plot(x,u_lf1(:,end)); hold on; plot(x,u_lw1(:,end));plot(x,u_up1(:,end)); plot(x, u_exact1(x,t(end)))
xlabel("x")
title("Plot of u(x,"+string(t(end))+")")
legend("Lax-Friedrichs","Lax-Wendroff", "Upwind", "Exact")


% For 'g_sq'

[~, ~, u_lf2] = hyperbolic1D(1, a, u_b2, N, l, D, T);

[~, ~, u_lw2] = hyperbolic1D(2, a, u_b2, N, l, D, T);

[~, ~, u_up2] = hyperbolic1D(3, a, u_b2, N, l, D, T);

figure
plot(x,u_lf2(:,end)); hold on; plot(x,u_lw2(:,end));plot(x,u_up2(:,end)); plot(x, u_exact2(x,t(end)))
xlabel("x")
title("Plot of u(x,"+string(t(end))+")")
legend("Lax-Friedrichs","Lax-Wendroff", "Upwind", "Exact")


%%
figure
plot1c1(0.05:0.05:0.5, u_b2, 500, 3);

figure
plot1c2(0.4, u_b2, 100:100:1500, 3);


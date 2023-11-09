%% Computer Exercise 4, David Ahnlund, Emil Gestsson
clc, clear variables;

global a D T

a = 2; D = 10; T = 4;
 
N = 200; 
l = 0.2;


tau = 2;
u_b1 = @(t) sin(2*pi*t/tau);

u_b2 = @(t) sign(sin(2*pi*t/tau));


% For 'g_sin'

[t, x, u_lf1] = hyperbolic1D(1, N, u_b1, l);

[~, ~, u_lw1] = hyperbolic1D(2, N, u_b1, l);

[~, ~, u_up1] = hyperbolic1D(3, N, u_b1, l);

plot(x,u_lf1(:,end)); hold on; plot(x,u_lw1(:,end));plot(x,u_up1(:,end))
xlabel("x")
legend("Lax-Friedrichs","Lax-Wendroff", "Upwind")

% For 'g_sq'

[~, x, u_lf2] = hyperbolic1D(1, N, u_b2, l);

[~, ~, u_lw2] = hyperbolic1D(2, N, u_b2, l);

[~, ~, u_up2] = hyperbolic1D(3, N, u_b2, l);

figure
plot(x,u_lf2(:,end)); hold on; plot(x,u_lw2(:,end));plot(x,u_up2(:,end))
xlabel("x")
legend("Lax-Friedrichs","Lax-Wendroff", "Upwind")


function [t, x, u] = hyperbolic1D(scheme, N, bound, lambda)

    global D T a dx dt t
    
    %scheme: 1 = Lax-Friedrichs, 2 = Lax-Wendroff, 3 = Upwind
    
    dx = D/N;
    dt = lambda*dx;
    t = 0:dt:T;
    x = 0:dx:D;

    u = zeros(length(x), length(t));
    u(1,:) = bound(t);
    
    
    if scheme == 1  % Lax-Friedrich
        k = dx^2;
    end
    if scheme == 2  % Lax-Wendroff
        k = a^2*dt^2;
    end
    if scheme == 3  % Upwind
        k = abs(a)*dx*dt;
    end

    for n = 2:length(t)-1
        for i = 2:length(x)-1

            u(i,n+1) = u(i,n) - a*lambda/2*(u(i+1,n)-u(i-1,n)) + 1/2*k*(u(i+1,n)-2*u(i,n)+u(i-1,n))/dx^2;

            if i == length(x)-1
                u(i+1,n) = 2*u(i,n)-u(i-1,n);
            end
        end
    end

end
%% Computer Exercise 4, David Ahnlund, Emil Gestsson
clc, clear variables;

global a D T

a = 2; D = 10; T = 4;
 
N = 100; 
l = 0.4;


tau = 2;
u_b1 = @(t) sin(2*pi*t/tau);

u_b2 = @(t) sign(sin(2*pi*t/tau));

[t,x,u] = hyperbolic1D(100,u_b2,l);

surf(t,x,u)

function [t, x, u] = hyperbolic1D(N, bound, l)

    global D T a dx dt t
    
    %scheme: 1 = Lax-Friedrichs, 2 = Upwind, 3 = Lax-Wendroff
    
    dx = D/N;
    dt = l*dx;
    t = 0:dt:T;
    x = 0:dx:D;

    u = zeros(length(x), length(t));
    disp(size(u))
    u(1,:) = bound(t);
    
    
    for n = 2:length(t)-1
        for i = 2:length(x)-1

            u(i,n+1) = u(i,n) - a*l/2*(u(i+1,n)-u(i-1,n)) + a^2*l^2/2*(u(i+1,n)-2*u(i,n)+u(i-1,n));

            if i == length(x)-1
                u(i+1,n) = 2*u(i,n)-u(i-1,n);
            end
        end
    end

end
%% Computer Exercise 4, David Ahnlund, Emil Gestsson
clc, clear variables;

global a D T

a = 2; D = 10; tau = 2; N = 100; T = 4;


u_b1 = @(t) sin(2*pi*t/tau);

u_b2 = @(t) sign(sin(2*pi*t/tau));

% For 'g_sin'

[t, x, u_lf1] = hyperbolic1D(1, 100, u_b1, 0.4);

[~, ~, u_up1] = hyperbolic1D(2, 100, u_b1, 0.4);

[~, ~, u_lw1] = hyperbolic1D(3, 100, u_b1, 0.4);

plot(x,u_lf1(:,end)); hold on; plot(x,u_up1(:,end));plot(x,u_lw1(:,end))
xlabel("x")
legend("Lax-Friedrich", "Upwind", "Lax-Wendroff")
% For 'g_sin'

[~, x, u_lf2] = hyperbolic1D(1, 100, u_b2, 0.4);

[~, ~, u_up2] = hyperbolic1D(2, 100, u_b2, 0.4);

[~, ~, u_lw2] = hyperbolic1D(3, 100, u_b2, 0.4);

figure
plot(x,u_lf2(:,end)); hold on; plot(x,u_up2(:,end));plot(x,u_lw2(:,end))
xlabel("x")
legend("Lax-Friedrich", "Upwind", "Lax-Wendroff")


function [t, x, u] = hyperbolic1D(scheme, N, bound, lambda)

    global D T a dx dt t
    
    %scheme: 1 = Lax-Friedrichs, 2 = Upwind, 3 = Lax-Wendroff
    
    dx = D/N;
    dt = lambda*dx;
    t = 0:dt:T;
    
    if scheme == 1
        A = 1/2 * ((1+lambda*a)*diag(ones(N-1,1),-1) + (1-lambda*a)*diag(ones(N-1,1),1));
    end
    if scheme == 2
        A = (a>0)*((a*lambda)*diag(ones(N-1,1),-1) + (1-a*lambda)*diag(ones(N,1),0))...
            +(a<0)*((-a*lambda)*diag(ones(N-1,1), 1) + (1+a*lambda)*diag(ones(N,1),0));
    end
    if scheme == 3
        A = (a^2*lambda^2/2+a*lambda/2)*diag(ones(N-1,1),-1) + (1-a^2*lambda^2)*diag(ones(N,1),0) + (a^2*lambda^2/2-a*lambda/2)*diag(ones(N-1,1),1);
    end
    
    u0 = zeros(N,1);
    u = zeros(length(u0), length(t));
    u(:,1) = u0;
    
    uk = u0;
    for n = 2:length(t)
        u_new = A*uk + [bound(t(n)); zeros(N-2,1); 2*uk(end)-uk(end-1)];
        u(:,n) = u_new;
        uk = u_new;
    end

    x = 0:dx:D-dx;

end
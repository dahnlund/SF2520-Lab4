%% Computer Exercise 4, David Ahnlund, Emil Gestsson
clc, clear variables;

a = 2;
D = 10;
tau = 2;
N = 100;
T = 4;

u_b1 = @(t) sin(2*pi*t/tau);

u_b2 = @(t) sign(sin(2*pi*t/tau));

%% Lax Friedrich

lambda = 0.4;
dx = D/N;
dt = lambda*dx;
t = 0:dt:T;
A = 1/2 * ((1+lambda*a)*diag(ones(N-1,1),-1) + (1-lambda*a)*diag(ones(N-1,1),1));

u0 = zeros(N,1);
saved_u = zeros(length(u0), length(t));
saved_u(:,1) = u0;

uk = u0;
for n = 2:length(t)
    u_new = A*uk + [u_b2(t(n)); zeros(N-2,1); 2*uk(end)-uk(end-1)];
    if ismember(n,2:1000)
    disp(2*uk(end)-uk(end-1))
    end
    saved_u(:,n) = u_new;
    uk = u_new;
end

x = 0:dx:D-dx;

surf(t,x,saved_u)
shading interp

%% Lax Wendroff

lambda = 0.4;
dx = D/N;
dt = lambda*dx;
t = 0:dt:T;
A = (a^2*lambda^2/2+a*lambda/2)*diag(ones(N-1,1),-1) + (1-a^2*lambda^2)*diag(ones(N,1),0) + (a^2*lambda^2/2-a*lambda/2)*diag(ones(N-1,1),1);

u0 = zeros(N,1);
saved_u = zeros(length(u0), length(t));
saved_u(:,1) = u0;

uk = u0;
for n = 2:length(t)
    u_new = A*uk + [u_b2(t(n)); zeros(N-2,1); 2*uk(end)-uk(end-1)];
    if ismember(n,2:1000)
    disp(2*uk(end)-uk(end-1))
    end
    saved_u(:,n) = u_new;
    uk = u_new;
end

x = 0:dx:D-dx;

surf(t,x,saved_u)
shading interp


surf(t,x,saved_u)
shading interp

a = 2; D = 10; T = 4;
 
N = 200; 
lambda = 0.2;
scheme = 1;
tau = 2;
u_b1 = @(t) sin(2*pi*t/tau);

bound = @(t) sign(sin(2*pi*t/tau));

%scheme: 1 = Lax-Friedrichs, 2 = Lax-Wendroff, 3 = Upwind

dx = D/N;
dt = lambda*dx;
t = 0:dt:T;
x = 0:dx:D;

u = zeros(length(x), length(t));
u(1,:) = bound(t);

if scheme == 1  % Lax-Friedrichs
    k = dx^2;
end
if scheme == 2  % Lax-Wendroff
    k = a^2*dt^2;
end
if scheme == 3  % Upwind
    k = abs(a)*dx*dt;
end

A = (a*lambda/2+k/(dx*2))*diag(ones(N-1,1),-1) + (1-k/dx)*diag(ones(N,1),0) + (a*lambda/2-k/(dx*2))*diag(ones(N-1,1),1);
A = [(a*lambda/2+k/(dx*2))*ones(N,1), A];
uk = u(:,1);
disp(uk)
for n = 1:length(t)-1
    u_new = A*uk;
    u(:,n) = u_new;
    uk = u_new;
end
u(:,end) = 2*u(end-1) - u(end-2);


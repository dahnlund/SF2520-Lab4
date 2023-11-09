%% SF2520 CE4, David Ahnlund, Emil Gestsson
clc, clear variables;

Fr = 0.35;
alpha = 1/Fr;
N = 100;
L0 = -0.4;
L1 = 0.7;
T = 0.15;

lambda = 0.1;

f = @(x,t) (sin(20*pi*x).*(abs(x) < 1/20))' * (sin(40*pi*t + pi/6)>0.5);


A= [1, alpha; alpha, 1];

F = @(x,t) [0;f(x,t)];
    
dx = abs(L0-L1)/N;
dt = lambda*dx;
t = 0:dt:T;
x = L0:dx:L1;

%% Upwind

[S, L] = eig(A); Lp = L.* (L>0); Lm = L.*(L<0);

Am = S*(S\Lm)'; Ap = S*(S\Lp)';

v = zeros(2,length(x), length(t));

for n = 2:length(t)-1
    for i = 2:length(x)-1
        
        if i == 2
            v(:,i-1,n) = 2*v(:,i,n)-v(:,i+1,n);
        end

        v(:,i,n+1) = v(:,i,n) - lambda*Ap*(v(:,i,n)-v(:,i-1,n)) - lambda*Am*(v(:,i+1,n)-v(:,i,n)) + F(x(i),t(n));
        
        if i == length(x)-1
            v(:,i+1,n) = 2*v(:,i,n)-v(:,i-1,n);
        end

    end
end

u1 = reshape(v(1,:,:), length(x), length(t));
u2 = reshape(v(2,:,:), length(x), length(t));
surf(t,x, u1)
figure
surf(t,x, u2)
figure

%% Laxen

v = zeros(2,length(x), length(t));

for n = 2:length(t)-1
    for i = 2:length(x)-1
        
        if i == 2
            v(:,i-1,n) = 2*v(:,i,n)-v(:,i+1,n);
        end
        
        F_tilde = (F(x(i), t(n)) + F(x(i), t(n+1)))/2 - lambda/4 * A * (F(x(i+1), t(n)) - F(x(i-1), t(n)));
        v(:,i,n+1) = v(:,i,n) - lambda/2*A*(v(:,i+1,n)-v(:,i-1,n)) + lambda^2/2*A^2*(v(:,i+1,n)-2*v(:,i,n)+v(:,i-1,n)) + F_tilde;
        
        if i == length(x)-1
            v(:,i+1,n) = 2*v(:,i,n)-v(:,i-1,n);
        end

    end
end

u1 = reshape(v(1,:,:), length(x), length(t));
u2 = reshape(v(2,:,:), length(x), length(t));
surf(t,x, u1)
figure
surf(t,x, u2)

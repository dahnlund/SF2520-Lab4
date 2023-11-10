%% SF2520 CE4, David Ahnlund, Emil Gestsson
function part2(CFL, N, plot_option, animation)

% CFL: ratio dt/dx, has to be < 1/(1+alpha), where alpha defined as below
% N: Discretization resoltion in x
% plot_option: if false, do not plot 3D
% animation: if true, animate the evolution of u in the Lax-Wendroff case

Fr = 0.35;
alpha = 1/Fr;
L0 = -0.4;
L1 = 0.7;
T = 0.15;

lambda = CFL;

f = @(x,t) (sin(20*pi*x).*(abs(x) < 1/20))' * (sin(40*pi*t + pi/6)>0.5);


A= [1, alpha; alpha, 1];

F = @(x,t) [0;f(x,t)];
    
dx = abs(L0-L1)/N;
dt = lambda*dx;
t = 0:dt:T;
x = L0:dx:L1;

%Upwind:

[S, L] = eig(A); Lp = L.* (L>0); Lm = L.*(L<0);

Am = S*(S\Lm)'; Ap = S*(S\Lp)';

v = zeros(2,length(x), length(t));

for n = 2:length(t)-1
    for i = 2:length(x)-1
        
        if i == 2
            v(:,i-1,n) = 2*v(:,i,n)-v(:,i+1,n);
        end

        v(:,i,n+1) = v(:,i,n) - lambda*Ap*(v(:,i,n)-v(:,i-1,n)) - lambda*Am*(v(:,i+1,n)-v(:,i,n)) + dt*F(x(i),t(n));
        
        if i == length(x)-1
            v(:,i+1,n) = 2*v(:,i,n)-v(:,i-1,n);
        end

    end
end

u_up = reshape(v(1,:,:), length(x), length(t));
v_up = reshape(v(2,:,:), length(x), length(t));

% Lax-Wendroff

v = zeros(2,length(x), length(t));

for n = 2:length(t)-1
    for i = 2:length(x)-1
        
        if i == 2
            v(:,i-1,n) = 2*v(:,i,n)-v(:,i+1,n);
        end
        
        F_tilde = (F(x(i), t(n)) + F(x(i), t(n+1)))/2 - lambda/4 * A * (F(x(i+1), t(n)) - F(x(i-1), t(n)));
        v(:,i,n+1) = v(:,i,n) - lambda/2*A*(v(:,i+1,n)-v(:,i-1,n)) + lambda^2/2*A^2*(v(:,i+1,n)-2*v(:,i,n)+v(:,i-1,n)) + dt*F_tilde;
        
        if i == length(x)-1
            v(:,i+1,n) = 2*v(:,i,n)-v(:,i-1,n);
        end

    end

    if animation == true
    u_lax = reshape(v(1,:,:), length(x), length(t));
    plot(x, u_lax(:,n), LineWidth=2)
    title("Lax-Wendroff, t = " + string(round(t(n),4)))
    ylabel("u")
    xlabel("x")
    xlim([L0 L1])
    ylim([-0.02 0.02])
    drawnow
    end
    
end

u_lax = reshape(v(1,:,:), length(x), length(t));
v_lax = reshape(v(2,:,:), length(x), length(t));

if plot_option ~= false
figure
mesh(t,x, u_lax)
title("Lax-Wendroff")
ylabel("x")
xlabel("t")

figure
mesh(t,x, u_up)
title("Upwind")
ylabel("x")
xlabel("t")

end

% Compare plots

figure
plot(x, u_up(:, end)); hold on; plot(x, u_lax(:,end));
title("u(x," + string(round(t(end),4))+")")
legend("Upwind", "Lax-Wendroff")
xlabel("x")
figure
plot(x, v_up(:,end)); hold on; plot(x, v_lax(:,end));
title("v(x," + string(round(t(end),4))+")")
legend("Upwind", "Lax-Wendroff")
xlabel("x")
end

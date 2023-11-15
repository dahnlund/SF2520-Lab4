using LinearAlgebra
using Plots

plotly()

CFL = 0.2
N = 300
Fr = 0.35
alpha = 1/Fr
L0 = -0.4
L1 = 0.7
T = 0.15

lambda = CFL;

function f(x,t)
    return (sin(20*pi*x).*(abs(x) < 1/20))' * (sin(40*pi*t + pi/6)>0.5)
end
function F(x,t)
    return [0;f(x,t)]
end
    
dx = abs(L0-L1)/N
dt = lambda*dx
t = 0:dt:T
x = L0:dx:L1

A = [1 alpha; alpha 1]

#Upwind:

L, S = eigen(A); L = diagm(L)

Lp = L.* (L.>0); Lm = L.*(L.<0)

Am = S*(S\Lm)'; Ap = S*(S\Lp)'

v = zeros(2,length(x), length(t));

for n in 2:length(t)-1
    for i in 2:length(x)-1
        
        if i == 2
            v[:,i-1,n] = 2*v[:,i,n]-v[:,i+1,n]
        end

        v[:,i,n+1] = v[:,i,n] - lambda*Ap*(v[:,i,n]-v[:,i-1,n]) - lambda*Am*(v[:,i+1,n]-v[:,i,n]) + dt*F(x[i],t[n])
        
        if i == length(x)-1
            v[:,i+1,n] = 2*v[:,i,n]-v[:,i-1,n]
        end

    end
end

u_up = reshape(v[1,:,:], length(x), length(t))
v_up = reshape(v[2,:,:], length(x), length(t))

# Lax-Wendroff

v = zeros(2,length(x), length(t));

for n = 2:length(t)-1
    for i = 2:length(x)-1
        
        if i == 2
            v[:,i-1,n] = 2*v[:,i,n]-v[:,i+1,n];
        end
        
        F_tilde = (F(x[i], t[n]) + F(x[i], t[n+1]))/2 - lambda/4 * A * (F(x[i+1], t[n]) - F(x[i-1], t[n]));
        v[:,i,n+1] = v[:,i,n] - lambda/2*A*(v[:,i+1,n]-v[:,i-1,n]) + lambda^2/2*A^2*(v[:,i+1,n]-2*v[:,i,n]+v[:,i-1,n]) + dt*F_tilde;
        
        if i == length(x)-1
            v[:,i+1,n] = 2*v[:,i,n]-v[:,i-1,n];
        end

    end  
end

u_lax = reshape(v[1,:,:], length(x), length(t));
v_lax = reshape(v[2,:,:], length(x), length(t));



plot1 = surface(t,x,u_up, size= (800,800), title = "Upwind", xlabel = "t", ylabel = "x")

plot2 = surface(t,x,u_lax, size = (800,800), title = "Lax Wendroff", xlabel = "t", ylabel = "x")
 
display(plot1)

display(plot2)

# Compare plots

plot_u_end = plot(x, u_up[:,end], label = "Upwind")
plot!(x, u_lax[:,end], label="Lax-Wendroff")
title!("u(x, " * string(round(t[end], digits  =4)) * ")")
xlabel!("x")

plot_v_end = plot(x, v_up[:,end],label="Upwind")
plot!(x, v_lax[:,end], label="Lax-Wendroff")
title!("v(x, " * string(round(t[end], digits = 4)) * ")")
xlabel!("x")

# Display the combined plots
display(plot_u_end)
display(plot_v_end)


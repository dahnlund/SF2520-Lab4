function hyperbolic(scheme, a, bound, N, lambda, D, T)
    
    #scheme: 1 = Lax-Friedrichs, 2 = Lax-Wendroff, 3 = Upwind

    dx = D/N
    dt = lambda*dx
    t = 0:dt:T
    x = 0:dx:D

    u = zeros(length(x), length(t))
    u[1,:] = bound(t)
    
    
    if scheme == 1  # Lax-Friedrichs
        k = dx^2
    end
    if scheme == 2  # Lax-Wendroff
        k = a^2*dt^2
    end
    if scheme == 3  # Upwind
        k = abs(a)*dx*dt
    end

    for n = 2:length(t)-1
        for i = 2:length(x)-1

            u[i,n+1] = u[i,n] - a*lambda/2*(u[i+1,n]-u[i-1,n]) + 1/2*k*(u[i+1,n]-2*u[i,n]+u[i-1,n])/dx^2

            if i == length(x)-1
                u[i+1,n] = 2*u[i,n]-u[i-1,n]
            end
        end
    end

    return t, x, u
end
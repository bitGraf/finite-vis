L = 1;
alpha = 111;

num_x = 100;
num_y = 100;

[x,y] = meshgrid(...
    linspace(0, L, num_x),...
    linspace(0, L, num_y));

N = 100; % num Y solns
M = 100; % num X solns
% Calc coefficients
A_mn = zeros(M, N);
lam2_mn = zeros(M, N);
mu_m = zeros(M,1);
mu_n = zeros(N,1);

t = 0;
dt = 0.000001;
while t < .001
    u = zeros(num_x,num_y);
    for m = 1:2:M
        mu_m(m) = -(pi/L)*m;
        for n = 1:N
            if mod(n,4) == 0
                continue;
            end
            mu_n(n) = -(pi/L)*n;

            x_eq = (-1)^m - 1;
            y_eq = cos(n*pi/2) - 1;
            A_mn(m,n) = 200/(m*n*pi^2) * x_eq*y_eq;

            lam2_mn(m,n) = alpha*(pi/L)^2*(m^2 + n^2);

            u = u + A_mn(m,n)*sin(mu_m(m)*x).*sin(mu_n(n)*y).*exp(-lam2_mn(m,n)*t);
        end
    end
    
    figure(1);
    clf;
    grid on;
    hold on;
    colormap(pink);
    view(-100, 22);
    xlabel('X (m)');
    ylabel('Y (m)');
    zlabel('t (K)');
    zlim([-10, 60]);
    s = surf(x,y,u);
    s.EdgeAlpha = 0.3;
    title(sprintf('T = %.3f ms', t*1000));
    colorbar;
    caxis([0,50]);
    drawnow
    %pause(.1);
%     break;

    t = t + dt;
    dt = dt*1.05;
end
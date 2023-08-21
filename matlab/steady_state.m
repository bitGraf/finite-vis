cmap = 'parula';

% define the domain
W = 1;
H = 1;

num_x = 200;
num_y = 300;

[x,y] = meshgrid(...
    linspace(0, W, num_x),...
    linspace(0, H, num_y));

% Boundary conditions
C1 = 100; % constant BC on one wall
C2 = 200;
C3 = 300;
C4 = 400;

% Coefficients
N = 100;

wn_x = @(n) -n*pi/W;
eta_y = @(n,y) exp(wn_x(n)*(2*H-y)) - exp(wn_x(n)*y);

wn_y = @(n) -n*pi/H;
eta_x = @(n,x) exp(wn_y(n)*(2*W-x)) - exp(wn_y(n)*x);

% solutions
u1 = zeros(num_y,num_x);
u2 = zeros(num_y,num_x);
u3 = zeros(num_y,num_x);
u4 = zeros(num_y,num_x);

for n = 1:1:N
    % u1(x,y) -> u1(x,0)=C1
    B_n = -2*C1*(cos(wn_x(n)*W) - 1) / (wn_x(n)*W*eta_y(n,0));
    u1 = u1 + B_n*eta_y(n,y).*sin(wn_x(n)*x);
    
    % u4(x,y) -> u4(0,y)=C4
    B_n = -2*C4*(cos(wn_y(n)*H) - 1) / (wn_y(n)*H*eta_x(n,0));
    u4 = u4 + B_n*eta_x(n,x).*sin(wn_y(n)*y);
    
    % u2(x,y) -> u2(W,y)=C2
    B_n = -2*C2*(cos(wn_y(n)*H) - 1) / (wn_y(n)*H*sinh(wn_x(n)*W));
    u2 = u2 + B_n*sinh(wn_x(n)*x).*sin(wn_y(n)*y);
    
    % u3(x,y) -> u2(x,H)=C3
    B_n = -2*C3*(cos(wn_x(n)*W) - 1) / (wn_x(n)*W*sinh(wn_y(n)*H));
    u3 = u3 + B_n*sinh(wn_y(n)*y).*sin(wn_x(n)*x);
end

us = u1 + u2 + u3 + u4;
tmax = max(max(us));
tmin = min(min(us));

figure(1);
clf;
grid on;
hold on;
axis square;
colormap(cmap);
shading interp;
view(-150, 22);
caxis([tmin, tmax]);
colorbar;
xlabel('X (m)');
ylabel('Y (m)');
zlabel('t (K)');
s = surf(x,y,us);
s.EdgeAlpha = 0;
s.FaceAlpha = 1;
plot3(xlim, 0*[1 1], C1*[1 1], 'k--', 'LineWidth', 4);
plot3(W*[1 1], ylim, C2*[1 1], 'k--', 'LineWidth', 4);
plot3(xlim, H*[1 1], C3*[1 1], 'k--', 'LineWidth', 4);
plot3(0*[1 1], ylim, C4*[1 1], 'k--', 'LineWidth', 4);
title('U_s(x,y) = U_1 + U_2 + U_3 + U_4');

%%
figure(2);
clf;
colormap(cmap);
shading interp;

subplot(2,2,1);
grid on;
hold on;
axis square;
view(-45, 33);
caxis([tmin, tmax]);
xlabel('X (m)');
ylabel('Y (m)');
zlabel('t (K)');
zlim([0,500]);
s = surf(x,y,u1);
s.EdgeAlpha = 0;
s.FaceAlpha = 1;
plot3(xlim, 0*[1 1], C1*[1 1], 'k--', 'LineWidth', 4);
title('U_1(x,y)');

subplot(2,2,2);
grid on;
hold on;
axis square;
view(-45, 33);
caxis([tmin, tmax]);
xlabel('X (m)');
ylabel('Y (m)');
zlabel('t (K)');
zlim([0,500]);
s = surf(x,y,u2);
s.EdgeAlpha = 0;
s.FaceAlpha = 1;
plot3(W*[1 1], ylim, C2*[1 1], 'k--', 'LineWidth', 4);
title('U_2(x,y)');

subplot(2,2,3);
grid on;
hold on;
axis square;
view(-45, 33);
caxis([tmin, tmax]);
xlabel('X (m)');
ylabel('Y (m)');
zlabel('t (K)');
zlim([0,500]);
s = surf(x,y,u3);
s.EdgeAlpha = 0;
s.FaceAlpha = 1;
plot3(xlim, H*[1 1], C3*[1 1], 'k--', 'LineWidth', 4);
title('U_3(x,y)');

subplot(2,2,4);
grid on;
hold on;
axis square;
view(-45, 33);
caxis([tmin, tmax]);
xlabel('X (m)');
ylabel('Y (m)');
zlabel('t (K)');
zlim([0,500]);
s = surf(x,y,u4);
s.EdgeAlpha = 0;
s.FaceAlpha = 1;
plot3(0*[1 1], ylim, C4*[1 1], 'k--', 'LineWidth', 4);
title('U_4(x,y)');

%%
figure(3);
clf;
grid on;
hold on;
colormap(parula(10));
caxis([0, 500]);
colorbar;
shading interp;
xlabel('X (m)');
ylabel('Y (m)');
contourf(x,y,us);
title('U_s(x,y) = U_1 + U_2 + U_3 + U_4');
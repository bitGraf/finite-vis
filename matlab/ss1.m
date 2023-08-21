W = 2;
H = 3;

num_x = 200;
num_y = 300;

[x,y] = meshgrid(...
    linspace(0, W, num_x),...
    linspace(0, H, num_y));

N = 100; % num
% Calc coefficients
C1 = 100; % constant BC on one wall

wn = @(n) -n*pi/W;
eta = @(n,y) exp(wn(n)*(2*H-y)) - exp(wn(n)*y);

u = zeros(num_y,num_x);
for n = 1:1:N
    B_n = -2*C1*(cos(wn(n)*W) - 1) / (wn(n)*W*eta(n,0));
    
    up = B_n*eta(n,y).*sin(wn(n)*x);
%     if any(any(up == Inf)) || any(any(isnan(up)))
%         n
%         break;
%     end
    u = u + up;
end

figure(1);
clf;
grid on;
hold on;
colormap('jet');
view(-150, 22);
xlabel('X (m)');
ylabel('Y (m)');
zlabel('t (K)');
s = surf(x,y,u);
s.EdgeAlpha = 0.1;
s.FaceAlpha = 1;
plot3(xlim, 0*[1 1], C1*[1 1], 'k--');

% xlim([0.5084    0.5794]);
% ylim([0,1]);
% zlim([94.6557  103.1719]);
% view(0,0);
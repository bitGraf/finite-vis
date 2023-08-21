% define the domain
W = 1;
N = 23;

k = 4;
num_x = k*1000;
x = linspace(0, k*W, num_x);

% Boundary conditions
C1 = 1; % constant BC on one wall

% Coefficients

wn_x = @(n) -n*pi/W;

% true function
s = mod(x,2)-1;
u = (s<0) - (s>0);

figure(1);
clf;
grid on;
hold on;
plot(x, .5*(u+1), 'k', 'LineWidth', 2);

b = zeros(1,N);

% fourier series
un = zeros(1, num_x);
for n = 1:1:N
    n;
    % u3(x,y) -> u2(x,H)=C3
    B_n = -2*C1*(W*cos(wn_x(n)*W) - 1) / wn_x(n);
%     B_n = -2*C1*(sin(wn_x(n)*W) - wn_x(n)*W*cos(wn_x(n)*W) - wn_x(n)*W) / wn_x(n)^2;
    b(n) = B_n;
    un = un + B_n*sin(wn_x(n)*x);
%     plot(x, B_n*sin(wn_x(n)*x));
end

plot(x, 0.5*(un+1), 'b');
title(sprintf('N = %.0f', N));
% plot(xlim, .85*[1 1], 'k--');
% plot(xlim, .15*[1 1], 'k--');

%%
% figure(2);
% plot(1:(N-1), b)
% points go from [0,N] i.e. N+1 points in each direction

N = 4;
r = 0.25;

g = -r;
f = 1 + 4*r;

U = zeros(N+1, N+1);
U(1,:)   = 100; % Top
U(:,N+1) = 200; % Right
U(N+1,:) = 300; % Bottom
U(:,1)   = 200; % Left

%%
% Identity
I = eye(N+1);

% diag([0, g, g, g, 0])
G = diag([0, repmat(g, 1, N-1), 0]);

% tridiagonal
F = eye(N+1);
for n = 2:(N)
    F(n,n-1) = g;
    F(n,n) = f;
    F(n,n+1) = g;
end

%%
F_inv = inv(F);

Cp     = zeros(N+1, N+1, N-2);
Bp     = zeros(N+1, N+1, N-1);
Bp_inv = zeros(N+1, N+1, N-1);
Dp     = zeros(N+1, N-1);
X      = zeros(N+1, N-1);

% precalcs:
Bp(:,:,1) = F;
Bp_inv(:,:,1) = inv(Bp(:,:,1));
Cp(:,:,1) = Bp_inv(:,:,1)*G;

for i = 2:(N-1)
    Bp(:,:,i) = F - G*Cp(:,:,i-1);
    Bp_inv(:,:,i) = inv(Bp(:,:,i));
    
    if i ~= (N-1)
        Cp(:,:,i) = Bp_inv(:,:,i)*G;
    end
end

%% Loop
U
for k = 1:2
    U_last = U;
    for i = 1:(N-1)
        R = U_last(2,:)';
        if i == 1
            D = R - G*(U_last(1,:)');
        elseif i == (N-1)
            D = R - G*(U_last(N+1,:)');
        else
            D = R;
        end

        if i == 1
            Dp(:,i) = Bp_inv(:,:,i)*D;
        else
            Dp(:,i) = Bp_inv(:,:,i)*(D - G*Dp(:,i-1));
        end
    end

    X(:,N-1) = Dp(:,N-1);
    U(N,2:N) = X(2:N,N-1)';

    for i = (N-2):-1:1
        X(:, i) = Dp(:,i) - Cp(:,:,i)*X(:,i+1);
        U(i+1,2:N) = X(2:N,i)';
    end

    U
end
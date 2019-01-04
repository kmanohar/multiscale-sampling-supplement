function tree = mrDMD(Xraw, dt, r, max_cyc, L,stackflag)
% function tree = mrDMD(Xraw, dt, r, max_cyc, L, stackflag)

% Inputs:
% Xraw      n by m matrix of raw data,
%           n measurements, m snapshots
% dt        time step of sampling
% r         rank of truncation
% max_cyc   to determine rho, the freq cutoff, compute
%           oscillations of max_cyc in the time window
% L         number of levels remaining in the recursion
%
% Code adapted from JN Kutz, SL Brunton, BW Brunton and JL Proctor,
% "Dynamic Mode Decomposition", SIAM
%
% Modified 2018/12/31

div = 2;

T = (size(Xraw, 2)) * dt; %compensate for matlab indexing

rho = max_cyc/T; % high freq cutoff at this level
sub = ceil(1/rho/8/pi/dt); % 4x Nyquist for rho
%sub = 1;
%% DMD at this level
Xaug = Xraw(:, 1:sub:end); % subsample
if (stackflag)
    Xaug = [Xaug(:, 1:end-1); Xaug(:, 2:end)];
end
X = Xaug(:, 1:end-1);
Xp = Xaug(:, 2:end);

% compute forward DMD
[U, S, V] = svd(X, 'econ');
if (r==0) % optimize projection rank at each level
    tf = isoutlier(diag(S));
    r = find(~tf,1);
else % user-specified rank 
    r = min(size(U,2),r);
end

U_r = U(:, 1:r); % rank truncation
S_r = S(1:r, 1:r);
V_r = V(:, 1:r);


Atilde = U_r' * Xp * V_r / S_r;

[W, D] = eig(Atilde);
lambda = diag(D);

Phi = Xp * V_r / S_r * (W/D);


%% compute power of modes
Vand = zeros(r, size(X, 2)); % Vandermonde matrix
for k = 1:size(X, 2),
    Vand(:, k) = lambda.^(k-1);
end;

% the next 5 lines follow Jovanovic et al, 2014 code:
G = S_r * V_r';
P = (W'*W).*conj(Vand*Vand');
q = conj(diag(Vand*G'*W));
Pl = chol(P,'lower');
b = (Pl')\(Pl\q); % Optimal vector of amplitudes b

%b = Phi\X(:,1);
%% consolidate slow modes, where abs(omega) < rho
omega = log(lambda)/sub/dt/2/pi;

mymodes = find(abs(omega) <= rho);

thislevel.T = T;
thislevel.rho = rho;
thislevel.hit = numel(mymodes) ;
thislevel.omega = omega(mymodes); 
thislevel.P = abs(b(mymodes));
thislevel.Phi = Phi(:, mymodes);
thislevel.lambda = lambda.^(1/sub); % properly scaled by sub

%% recurse on partitions
nextlevel = cell(div,1);
if L > 1,
    sep = floor(size(Xraw,2)/div);
    start = 1;
    for d=1:div
       nextlevel{d} = mrDMD(Xraw(:,start:start+sep-1),dt,r,max_cyc,L-1,stackflag);
       start = start+sep;
    end
%     nextlevel1 = mrDMD(Xraw(:, 1:sep),dt, r, max_cyc, L-1,stackflag);
%     nextlevel2 = mrDMD(Xraw(:, sep+1:end),dt, r, max_cyc, L-1,stackflag);
else
    for d=1:div
        nextlevel{d} = cell(0);
    end
%     nextlevel1 = cell(0);
%     nextlevel2 = cell(0);
end;

%% reconcile indexing on output
% (because Matlab does not support recursive data structures)
tree = cell(L, div^(L-1));
tree{1,1} = thislevel;

for l = 2:L,
    col = 1;
    for d=1:div        
        nxt = nextlevel{d};
        for j = 1:div^(l-2),
            tree{l, col} = nxt{l-1, j};
            col = col + 1;
        end
    end;
        

%     for j = 1:2^(l-2),
%        tree{l, col} = nextlevel1{l-1, j};
%        col = col + 1;
%     end;
%     for j = 1:2^(l-2)
%        tree{l, col} = nextlevel2{l-1, j};
%        col = col + 1;
%     end;
end;
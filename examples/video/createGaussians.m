function [ X, dt, modes ] = createGaussians( dist )
%createGaussians construct data for multiscale video example
%  Inputs:
%   dist	distance from the frid center of the two smaller modes
%
%  Outputs:
%   X       output matrix where each column is one reshaped time snapshot
%   dt      timestep 
%   modes	true spatial modes of the data
%  See script video_mrdmd_demo.m for usage
% 
% Code adapted from JN Kutz, SL Brunton, BW Brunton and JL Proctor,
% "Dynamic Mode Decomposition", SIAM
% Modified 2018/12/31

% DATASET PARAMETERS
nx = 80; % horizontal pixels
ny = 80; % vertical pixels
n  = nx * ny;
T = 10; % Total time duration in seconds
dt = 0.01; % Time step
t = dt:dt:T;

% 3 underlying modes
[Xgrid, Ygrid] = meshgrid(1:nx, 1:ny);

% center Gaussian 
mode1.u = exp(-((Xgrid-40).^2/250+(Ygrid-40).^2/250));
mode1.f = 5.55; % cycles/sec
mode1.A = 1;
mode1.lambda = exp(1i * mode1.f*2*pi*dt);
mode1.range = [0, 5];

leftcoord = 40-dist;
mode2.u = exp(-((Xgrid-leftcoord).^2/50+(Ygrid-leftcoord).^2/50));
mode2.f = 0.9; % cycles/sec 
mode2.A = 2;
mode2.lambda = exp(1i * mode2.f*2*pi*dt);
mode2.range = [2.5,7.5];

rightcoord = 40+dist;
mode3.u = exp(-((Xgrid-rightcoord).^2/50+(Ygrid-rightcoord).^2/50));
mode3.f = .15; % cycles/sec 
mode3.A = 1;
mode3.lambda = exp(1i * mode3.f*2*pi*dt);
mode3.range = [0, T];

modes(1) = mode1;
modes(2) = mode2;
modes(3) = mode3;

% make the movie
Xclean = zeros(n, numel(T));
Qoi = zeros(size(Xclean));

Snap = zeros(nx*ny, 3); % reusable temp array with 3 modes separated

tseries = zeros(3,length(t));

for ti = 1:numel(t),
    for mi = 1:numel(modes),        
        mymode = modes(mi);
        if ti > round(mymode.range(1)/dt) && ...
            ti < round(mymode.range(2)/dt),
            Snap(:,mi) = reshape(mymode.A * mymode.u * ...
                real(mymode.lambda^ti),nx*ny,1);
            tseries(mi,ti) = mymode.A*real(mymode.lambda^ti);
        end;
    end;
    Qoi(:, ti) = Snap(:,2);
    Xclean(:, ti) = sum(Snap,2);
end;

X = Xclean;

%% UNCOMMENT TO OUTPUT TIME SERIES OF DATA
% for i=1:3
%     plot(t,tseries(i,:),'k-','LineWidth',1.25);
%     set(gcf,'Position',[380 320 540 200]);
%     export_fig(['../../figures/video_tseries_' num2str(i)],'-pdf');
% end

end


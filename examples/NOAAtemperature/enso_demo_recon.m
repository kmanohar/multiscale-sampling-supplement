%% Script to reconstruct snapshots of NOAA SST from optimal multiscale sensors
% See Manohar, Kaiser, Brunton and Kutz, "Optimized Sampling 
% for Multiscale Dynamics", https://arxiv.org/abs/1712.05085
%
% Set PRINT_FIG=true to export figures into figpath 
% Output figures are saved into figpath
%
% Modified 2018-12-31

clear; close all; clc
datpath = 'data/';
figpath = '../../figures/';
f_input = 'enso';

% Set PRINT_FIG=true to export figures
PRINT_FIG = false;
set(0, 'defaultfigurecolor', 'w');

% Sea surface temperature NOAA
[Lat,Lon,time,mask,dat] = read_data_enso([datpath 'sst.wkmean.1990-present_2.nc'],...
    [datpath 'lsmask.nc']);


% 16 year period starting exactly beginning of 2001
% timeval = time(end-8-835+831+1:end);
% time = time(end-8-835:end-8-835+831);
% dat = dat(:,:,end-8-835:end);

numyears = 16;
nweeks = numyears*52;


Y = zeros(length(mask(mask==1)),size(dat,3));
for i=1:size(dat,3)
    sst = dat(:,:,i);
    Y(:,i) = sst(mask==1);           
end


[N,M] = size(Y);

% mrDMD frequencies should be yearly
% scale omegas by 52
% but keep dt=1 week.
dt = mean(diff(time))/7; 
timeval = time(nweeks+1:end);

% 16 year period ending exactly in 2017
time2 = time(end-8-831:end-8);
tree2 = mrDMD_fb(Y(:,end-8-831:end-8),dt,10,1,4, true);

% mrDMD, 16 year period beginning in 1990
time = time(1:nweeks);
tree = mrDMD_fb(Y(:,1:length(time)),dt,10,1,4, true);

[U,S,V] = svd(Y(:,1:length(time)),'econ');
[~,~,piv] = qr(U(:,1:30)',0);
qdeim = piv(1:30);

%% plot amplitudes of 1990 mrDMD

indt = [2  6  9 12]*52+26; % measure in summer months

[ptree, map, low_f_cutoff] = mrDMD_map(tree);
[L, J] = size(tree);

imagesc(-sqrt(map));

shortInt = tree{L,J}.T;

T = datetime(1800,1,1,0,0,0) + days(time(1:shortInt:end));
set(gca, 'YTick', 0.5:(L+0.5), 'YTickLabel', low_f_cutoff*52);
set(gca,'XTick', (0:J-1) + 0.5,'XTickLabel',num2str([month(T),year(T)],'%d/%d'));
axis xy; axis tight; box on; grid on

ylim = get(gca,'ylim');

times = year(datetime(1800,1,1,0,0,0) + days(time(indt)));
hold on
for i=1:length(indt)
    stem(1+indt(i)/shortInt, ylim(2),'k.','LineWidth',1.5)
    
end

set(gca, 'XTickLabelRotation',45,'ylim',ylim,'FontSize',12);
colormap(gca,'pink');shading flat

if PRINT_FIG
    export_fig([figpath 'mrdmd_recon_times'],'-png');
end

%% collect unique mrDMD modes

lmap = []; jmap = [];
Lib = [];
Omega = [];
Amp = [];
Periods = [];

for l=1:L
    for j=1:2^(l-1)
        Lib = [Lib tree{l,j}.Phi(1:N,:)];
        Omega = [Omega; tree{l,j}.omega]; % weekly
        Amp = [Amp; tree{l,j}.P];
        Periods = [Periods; tree{l,j}.T];
    end
end

[~,~,piv] = qr(Lib.','vector');
sens = piv(1:size(Lib,2));

Phi = cell(1,L);
omega = cell(L,1);
b = cell(L,1);
t0 = cell(L,1);

PhiJ = cell(1,J); win =1;
lab = [];
xL = zeros(N,L); %prediction at each level
update = false;
Lev_coeff = [];

Xrecon = [];

%% compute reconstruction of SST from multiscale sensor measurements
for i=1:length(time)
    x = Y(:,i);
    
    %xavg = L(:,i); % filtered data for sensors
    
    if (i>2 && i < length(time)-1)
        xavg = mean(Y(:,i-2:i+2),2);
    else 
        xavg = x;
    end


    % at each level update t-t0
    for l=1:L
        Int = tree{l,1}.T;
        
        if (mod(i,Int)==0 && i< length(time))
            t0{l} = time(i+1);
            ind = floor((i+1)/Int)+1;
            update = true;
        elseif (i==1)
            ind = 1;
            t0{l} = time(1);
            update = true;
        end
        
        if (update)   
            %disp([i ind]) % update modes when time windows change

            Phi{l} = tree{l,ind}.Phi(1:N,:);
            omega{l} = tree{l,ind}.omega;
            b{l} = tree{l,ind}.P;
            %sub{l} = ceil(1/8/pi/tree{l,ind}.rho/dt);
            

            Modes = cell2mat(Phi);
            
            ti = (time(i)-t0{l})/7; % time in weeks
            xL(:,l) = Phi{l}*(exp(2*pi*omega{l}*ti).*b{l});
            if (i==1)
                x0 = x;
            else
                x0 = Y(:,i+1);
            end
            lev_coeff = xL\x0;           
            
            if (l==L) % store J sets of modes by last window
                PhiJ{win} = cell2mat(Phi);
                lab = [lab; win*ones(size(PhiJ{win},2),1)];
                Lev_coeff = [Lev_coeff lev_coeff];
                win = win+1;
            end
            
            update = false;
            r = size(Modes,2);              
            disp([size(Modes,2) length(sens)]);
        end                
            
        ti = (time(i)-t0{l})/7; % time in weeks
        xL(:,l) = Phi{l}*(exp(2*pi*omega{l}*ti).*b{l});
    end
    
    
    % weighted mrDMD reconstruction
    % least-squares fit of by-level approximation to true snapshot
    xpred = xL*lev_coeff;    
    
    % collect Phi at all levels and approximate
    xhat = Modes*(Modes(sens,:)\x(sens));
    xhat2 = Lib*(Lib(sens,:)\x(sens));  
    
    xpod = U(:,1:r)*(U(qdeim,1:r)\x(qdeim));

    if (ismember(i, indt))
        Xrecon = [Xrecon xhat];
    end
    errpred(i) = norm(x-xpred)/norm(x);
    errsens(i) = norm(x-xhat)/norm(x);
    errsens2(i) = norm(x-xhat2)/norm(x);
    errpod(i) = norm(x-xpod)/norm(x);
    
end



%% display reconstructions of sea surface temperature from multiscale sensors

[~,~,time,mask,~] = read_data_enso([datpath 'sst.wkmean.1990-present_2.nc'],...
    [datpath 'lsmask.nc']);
times = datetime(1800,1,1,0,0,0) + days(time(indt));

for i=1:size(Xrecon,2)
    figure; 
    t = times(i);
    display_fig_sst(Xrecon(:,i),mask,[],[-1.8 35.6]);
    if PRINT_FIG
        export_fig([figpath 'recon' num2str([month(t),year(t)],'%d_%d')],'-pdf');
    end
    title(num2str([month(t),year(t)],'%d/%d'));
end

%% display sensors
figure;

display_sensors(mask,sens,'r');
if PRINT_FIG
    export_fig([figpath 'recon_sensors.pdf']);
end
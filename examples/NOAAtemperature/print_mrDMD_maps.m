%% Script to display mrDMD mode amplitude maps of NOAA sea surface temperature
% data, as shown in Manohar, Kaiser, Brunton and Kutz, "Optimized Sampling 
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
set(0,'defaultaxescolororder',parula(5));

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

% mrDMD, 16 year period ending exactly in 2017
time2 = time(end-8-831:end-8);
tree2 = mrDMD_fb(Y(:,end-8-831:end-8),dt,10,1,4, true);

% mrDMD, 16 year period beginning in 1990
time = time(1:nweeks);
tree = mrDMD_fb(Y(:,1:length(time)),dt,10,1,4, true);


%% plot amplitudes of 1990 mrDMD 
figure;
[ptree, map, low_f_cutoff] = mrDMD_map(tree);
[L, J] = size(tree);

imagesc(-sqrt(map));

shortInt = tree{L,J}.T;
T = datetime(1800,1,1,0,0,0) + days(time(1:shortInt:end));
set(gca, 'YTick', 0.5:(L+0.5), 'YTickLabel', low_f_cutoff*52);
set(gca,'XTick', (0:J-1) + 0.5);
axis xy; axis tight; box on; grid on
set(gca,'XTickLabel', strcat(num2str(month(T)),'/', num2str(year(T))));
set(gca, 'XTickLabelRotation',45);

colormap(gca,'pink');shading flat

if PRINT_FIG
    export_fig([figpath 'mrdmd_enso_2016_L4'],'-pdf');
	export_fig([figpath 'mrdmd_enso_1997'],'-png');
end

%% plot amplitudes of 2016 mrDMD 

[ptree, map, low_f_cutoff] = mrDMD_map(tree2);
[L, J] = size(tree);

imagesc(-sqrt(map));

shortInt = tree2{L,J}.T;
T = datetime(1800,1,1,0,0,0) + days(time2(1:shortInt:end));
set(gca, 'YTick', 0.5:(L+0.5), 'YTickLabel', low_f_cutoff*52);
set(gca,'XTick', (0:J-1) + 0.5);
axis xy; axis tight; box on; grid on
set(gca,'XTickLabel', strcat(num2str(month(T)),'/', num2str(year(T))));
set(gca, 'XTickLabelRotation',45);

colormap(gca,'pink');shading flat

if PRINT_FIG
    export_fig([figpath 'mrdmd_enso_2016'],'-png');
end

%% collect and display unique mrDMD modes

lmap = []; jmap = [];
Lib = []; Lib2 = [];
Omega = []; Omega2 = [];
Amp = []; Amp2 = [];
Periods = [];
Lambda = []; Lambda2 = [];

tol = 1e-2;

for l=1:L
    for j=1:2^(l-1)
        Lib = [Lib tree{l,j}.Phi(1:N,:)];        
        Omega = [Omega; tree{l,j}.omega*52]; % yearly
        Amp = [Amp; tree{l,j}.P];
        Periods = [Periods; tree{l,j}.T];
        Lambda = [Lambda; tree{l,j}.lambda];
        
        ind = find(imag(tree{l,j}.omega*52)>tol);
        for kk=1:length(ind)
            figure;
            display_fig_sst(abs(tree{l,j}.Phi(1:N,ind(kk))),mask,[],[]);
            str = ['Phi97' num2str(l) ',' num2str(j) '_' ...
                num2str(round(imag(tree{l,j}.omega(ind(kk))*52),2))];
            if PRINT_FIG
                export_fig([figpath str],'-pdf');
            end
        end
            
        Lib2 = [Lib2 tree2{l,j}.Phi(1:N,:)];
        Omega2 = [Omega2; tree2{l,j}.omega*52]; % yearly
        Amp2 = [Amp2; tree2{l,j}.P];
        Lambda2 = [Lambda2; tree2{l,j}.lambda];
        
        ind = find(imag(tree2{l,j}.omega*52)>tol);
        for kk=1:length(ind)
            figure;
            display_fig_sst(abs(tree2{l,j}.Phi(1:N,ind(kk))),mask,[],[]);
            str = ['Phi16' num2str(l) ',' num2str(j) '_' ...
                num2str(round(imag(tree2{l,j}.omega(ind(kk))*52),2))];
            if PRINT_FIG
                export_fig([figpath str],'-pdf');
            end
        end
    end
end

%% Filter library to keep only oscillatory and background mode
ind = find(abs(imag(Omega))>tol);
%ind = [1 ;ind];

ind2 = find(abs(imag(Omega2))>tol);
%ind2 = [1; ind2];

%plot(Omega(ind),'bo');hold on,plot(Omega2(ind2),'ro')

figure;
hold on
idx = find(abs(imag(Omega))>.1);
scatter(real(Omega),imag(Omega),150,0.5*[1 1 1])
scatter(real(Omega(idx)),imag(Omega(idx)),150,'r','filled');

grid on; box on

if PRINT_FIG
    export_fig([figpath 'mssp_spectrum.pdf']);
end



%% display QDEIM sensor locations and corresponding time series

r = 3;
cmap = lines(r);
[Phi,S,V] = svd(Y(:,1:length(time)),'econ');
Phi = Phi(:,1:r);
[~,~,piv] = qr(Phi','vector');
piv = piv(1:r);

figure;
display_sensors(mask,piv,cmap);box on
if PRINT_FIG
    export_fig([figpath 'qdeim_sensors.pdf']);
end

figure;
plot(Y(piv,1:length(time))','LineWidth',2); axis tight
shortInt = tree{L,J}.T;
T = datetime(1800,1,1,0,0,0) + days(time(1:shortInt:end));
set(gca,'XTick', (0:J-1)*shortInt + 0.5);
axis tight; box on;
set(gca,'XTickLabel', num2str(year(T)));
%set(gca,'XTick',[]);set(gcf,'Position',[380   320   540   200]);

if PRINT_FIG
    export_fig([figpath 'qdeim_tseries.pdf']);
end

%% display multiscale sensor locations and corresponding time series
cmap = lines(r);

Phi = Lib(1:N,idx);
[~,~,piv] = qr(Phi','vector');
piv = piv(1:r);

figure;
display_sensors(mask,piv,cmap); box on
if PRINT_FIG
    export_fig([figpath 'mssp_sensors.pdf']);
end


figure;
plot(Y(piv,1:length(time))','LineWidth',2); axis tight; box on

set(gca,'XTick', (0:J-1)*shortInt + 0.5);
axis tight; box on;
set(gca,'XTickLabel', num2str(year(T)));

%set(gca,'XTick',[]);set(gcf,'Position',[380   320   540   200]);

if PRINT_FIG
    export_fig([figpath 'mssp_tseries.pdf']);
end

%% display POD energy spectrum 

sig = diag(S)/sum(diag(S));

figure;
hold on
semilogy(sig,'.','Color',.5*[1 1 1]); 
scatter(1:r,sig(1:r),150,'r','filled');

set(gca,'yscale','log'); axis tight
grid on, box on
if PRINT_FIG
    export_fig([figpath 'POD_spectrum.pdf']);
end



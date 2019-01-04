%% Script for estimating El Nino coefficients using NOAA SST data
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

% Set PRINT_FIG=true to export figures
PRINT_FIG = false;
% Sea surface temperature NOAA
[Lat,Lon,time,mask,dat] = read_data_enso([datpath 'sst.wkmean.1990-present_2.nc'],...
    [datpath 'lsmask.nc']);

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


% mrDMD, 16 year period beginning in 1990
time = time(1:nweeks);
meanY = mean(Y(:,1:length(time)),2);
tree = mrDMD(Y(:,1:length(time)),dt,10,1,4, true);


%% plot amplitudes of 1990 mrDMD
figure;
[ptree, map, low_f_cutoff,Phi1] = mrDMD_map(tree);
[L, J] = size(tree);

imagesc(-sqrt(map));

shortInt = tree{L,J}.T; % each year
T = datetime(1800,1,1,0,0,0) + days(time(1:shortInt:end));
set(gca, 'YTick', 0.5:(L+0.5), 'YTickLabel', low_f_cutoff*52);
set(gca,'XTick', (0:J-1) + 0.5);
axis xy; axis tight; box on
set(gca,'XTickLabel', strcat(num2str(month(T)),'/', num2str(year(T))));
set(gca, 'XTickLabelRotation',45);
grid on

colormap(gca,'pink'); shading flat

if PRINT_FIG
    export_fig([figpath 'mrdmd_enso_1997'],'-png');
end

%% generate matrix library from mrDMD tree of modes

lmap = []; jmap = [];
Lib = []; 
Omega = []; 

for l=1:L
    for j=1:2^(l-1)
        Lib = [Lib tree{l,j}.Phi(1:N,:)];        
        Omega = [Omega; tree{l,j}.omega*52]; % yearly
    end
end



%% choose sensors only from identified slow modes
Phi1 = Phi1(1:N,:);
for k = 1:size(Phi1,2)
    Phi1(:,k) = Phi1(:,k)/norm(Phi1(:,k),2);
end

[~,~,pivk] = qr(Phi1.',0);
sens = pivk(1:size(Phi1,2));

display_sensors(mask,sens,'r');
%export_fig([figpath 'enso_demo_1_sensors'], '-pdf');

%% compute sparse library coefficients from optimal sensor measurements

A = [];
X = Y;

J = 832;
opts = spgSetParms('verbosity',0);         % Turn off the SPGL1 log output

for i=1:size(X,2)
    a = spg_bp(real(Lib(sens,:)),X(sens,i),opts);
    A = [A a];
end

shortInt = 52;

%% envelope of ENSO coefficient from training window 1990-2006

% Reload full temporal vector
[Lat,Lon,time,mask,dat] = read_data_enso([datpath 'sst.wkmean.1990-present_2.nc'],...
    [datpath 'lsmask.nc']);

% figure formatting options
set(0, 'defaultfigureposition', [380   320   540   200]);
set(0, 'defaultfigurecolor', 'w');

figure; 
plot(real((A(19,:))),'LineWidth',1.5)
[yu,yl] = envelope(real(A(19,:)),40,'peak');
ntt = length(yu);
axis tight
hold on, plot(1:ntt,yu,'k-',1:ntt,yl,'k-','LineWidth',.5,'Color',0.5*[1 1 1]);
%hold on, plot( mean([yu ;yl]),'LineWidth',1)
envavg = mean([yu;yl]); baseval = mean([yu yl]);
set(gca,'XTick',(1:52:length(time)-1));

times = datetime(1800,1,1,0,0,0) + days(time(2:shortInt:end)+1);
set(gca,'XTickLabel', strcat(num2str(month(times)),'/', num2str(year(times))));
set(gca,'XTickLabelRotation',45)
axis tight
set(gca,'xlim',[52,1390])

if PRINT_FIG
    export_fig([figpath 'enso_demo_1_env'],'-pdf')
end

%% prediction envelope mean

figure;
nino = envavg; nino(nino< baseval) = baseval;
nina = envavg; nina(nina> baseval) = baseval;
h = area(nino, 'BaseValue',baseval);
h(1).FaceColor = 'r'; hold on,
h = area(nina, 'BaseValue',baseval);
h(1).FaceColor = 'b'; 
grid on

set(gca,'XTick',(1:52:length(time)-1));

times = datetime(1800,1,1,0,0,0) + days(time(2:shortInt:end)+1);
set(gca,'XTickLabel', strcat(num2str(month(times)),'/', num2str(year(times))));
set(gca,'XTickLabelRotation',45)
axis tight
set(gca,'xlim',[52,1390])

if PRINT_FIG
    export_fig([figpath 'enso_demo_1_baseline'],'-pdf')
end

%% actual ENSO mode coefficient running means

a = A(19,:);
windowlen = 13; %3 months ~ 13 weeks
shift = 1;

runmean = zeros(length(a)-windowlen,1);
tt = zeros(size(runmean));
ntt = length(tt);

for i=1:(length(A)-windowlen)
    runmean(i) = mean(real(a(i:i+windowlen-1)));
    %each mean has a midpoint in time associated with it
    tt(i) = time(median(i:i+windowlen-1));
end

figure;% plot(runmean,'LineWidth',1.5);ba
[yu,yl] = envelope(runmean,40,'peak');
axis tight;

plot(1:ntt,runmean, '-','LineWidth',1.5);
axis tight
hold on, plot(1:ntt,yu,'k-',1:ntt,yl,'k-','LineWidth',.5,'Color',0.5*[1 1 1]);
hold on, plot(1:ntt, mean([yu yl],2),'LineWidth',1)
vec = mean([yu yl],2);
plot(get(gca,'xlim'), mean(vec)*[1 1],'r-')
set(gca,'XTick',(0:52:2*J-1) + 0.5);

times = datetime(1800,1,1,0,0,0) + days(tt(1:shortInt:end));
set(gca,'XTickLabel', strcat(num2str(month(times)),'/', num2str(year(times))));
set(gca,'XTickLabelRotation',45)
axis tight
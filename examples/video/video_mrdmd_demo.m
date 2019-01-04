%% Script to estimate temporal coefficients of mrDMD modes for multiscale video data
% Set PRINT_FIG=true to export figures into figpath 
% Output figures are saved into figpath
%
% Modified 2018-12-31

clear; close; clc
PRINT_FIG = false;
figpath = '../../figures/';

%% GENERATE DATASET 
% data defined on nx x ny grid
nx = 80;
ny = 80;

xc = 10; yc = 10;
[Xgrid, Ygrid] = meshgrid(1:nx, 1:ny);

F = exp(-.05*(Xgrid-41).^2-.05*(Ygrid-41).^2);
sig = [0 logspace(-10,0,11)];

[Xclean,dt,modes] = createGaussians(30); % Centers of Gaussian spatial modes (10,10), (40,40), (70,70)

%% UNCOMMENT TO GENERATE VIDEO FILE OF DATA
% gmin = min(Xclean(:));gmax = max(Xclean(:));
% vidObj = VideoWriter('video','MPEG-4');
% vidObj.FrameRate=60;
% open(vidObj);
% for i=1:size(Xclean,2)
%     frame = Xclean(:,i);    
%     frame = reshape(frame,nx,ny);
%     imagesc(frame,[gmin gmax]); colormap bone; 
%     axis xy, axis square, 
%     set(gca,'XTick',[],'YTick',[])
%     %set(gca,'nextplot','replacechildren');
%     xlabel(num2str(dt*i));    
%     currFrame = getframe;
%     writeVideo(vidObj,currFrame);
% end
%   
% % Close the file.
% close(vidObj);

%% COMPUTE MULTIRESOLUTION DMD FROM DATA
Phi = zeros(nx*ny,length(modes));

for i=1:length(modes)
    Phi(:,i) = reshape(modes(i).u,nx*ny,1);    
end

% mrDMD 
L = 6; % number of levels
r = 5; % rank of truncation
N = size(Xclean,1);
mrdmd = mrDMD(Xclean, dt, r, 2, L,true);

% compile visualization of multiresolution mode amplitudes
[ptree, map, low_f] = mrDMD_map(mrdmd);
[L, J] = size(mrdmd); T = 10;

Mr = id_mrDMD_modes( mrdmd, [.15;.9; 5.55], ptree);

% order by energy content
Mr = Mr(1:N,[3 2 1]);

%% COMPUTE MULTISCALE SENSOR PLACEMENTS
tic;
[~,~,pivot] = qr(Mr.',0);
sensors = pivot(1:3);
toc

% POD-DEIM sensor placement
[U,~,~] = svd(Xclean,'econ');
Ur = U(:,1:3);
[~,~,piv]= qr(Ur',0);
qdeim = piv(1:3);   
    
runs = 50;
%% COMPUTE TEMPORAL COEFFICIENTS FOR MRDMD & POD MODES FROM SENSORS

rng(971);

Epod = zeros(runs,length(sig));
E = zeros(runs,length(sig));

% Generate temporal coefficients using random process
Coefs = randn(3,runs);
for i=1:3
    Coefs(i,:) = conv(Coefs(i,:),ones(1,10)/15,'same'); %// apply a simple low-pass filter
    % normalize all coefficient rows
    Coefs(i,:) = Coefs(i,:)/sqrt(sum(Coefs(i,:).^2));
end

for i = 1:length(sig)
              
    A = zeros(3,runs);
    Apod = zeros(3,runs);
    
    for cv = 1:runs
        coefs = Coefs(:,cv);
        
        Phin = Phi*coefs;
        Phin = Phin + randn(size(Phin))*sig(i);

        a = Mr(sensors,:)\Phin(sensors,:); 
        A(:,cv)=a;
        apod = Ur(qdeim,:)\Phin(qdeim,:); 
        Apod(:,cv)=apod;
    end    
    
    % normalize for visual comparison    
    for k=1:3  
        [~,ind] = max(abs(Coefs(k,:)));
        val = Coefs(k,ind);
        A(k,:) = A(k,:)/sqrt(sum(A(k,:).^2));
        A(k,:) = sign(val*A(k,ind))*A(k,:);
        Apod(k,:) = Apod(k,:)/sqrt(sum(Apod(k,:).^2));
        Apod(k,:) = sign(val*Apod(k,ind))*Apod(k,:);
        
        if (i==9)
            figure;
            plot(Coefs(k,:),'ks-','MarkerSize',10), hold on
            plot(A(k,:),'.:'); plot(Apod(k,:),'.:'); axis tight
            if (k==3)
                legend('True','mrDMD','POD','Location','northeast');
            end
            
            set(gcf,'Position',[380 320 540 200]);
            
            %export_fig([figpath 'video_time_coefs' num2str(k)],'-pdf','-transparent');
        end
    end
    
    % 1-norm error metric indicative of visual comparison
    E(:,i) = sqrt(sum((Coefs-A).^2,1))';
    Epod(:,i) = sqrt(sum((Coefs-Apod).^2,1))';
    
    disp([median(Epod(:,i)) median(E(:,i))])
end

%save('video_noise_study','sig','Epod','E');

%% PLOT BOXPLOT OF TEMPORAL COEFFICIENT ERROR

figure;

cmap = get(0,'defaultAxesColorOrder');
% boxplot of reconstruction error 
dmdpod = [repmat({'  '},1,length(sig)),repmat({' '},1,length(sig))];
colorgrp = repmat({'1','2'},1,length(sig));
boxplot(abs([Epod E]),{repmat(sig,1,2),dmdpod},'Symbol','k.', ...
    'ColorGroup',colorgrp,'factorgap',[5 2],'labelverbosity','minor',...
    'Colors',cmap([2 1],:))

set(findobj(gca,'Type','text'),'FontSize',14,'VerticalAlignment','top',...
    'Rotation',45); 
set(findobj(gca,'Type','line'),'linew',1);
%set(findobj(gca,'Type','line','-and','Tag','Box','-and','Color',[.75 0 0]),'Color',cmap(1,:));
%set(findobj(gca,'Type','line','-and','Tag','Box','-and','Color',[0 .75 .75]),'Color',cmap(2,:));

set(gcf,'Position',[380 320 600 250]);
set(gca,'FontSize',14,'yscale','log');
legend(findall(gca,'Tag','Box'), {'mrDMD','POD'},'Location','best');

if PRINT_FIG
    export_fig([figpath 'video_noise_study'], '-pdf');
end

%% PLOT MRDMD MODE AMPLITUDE MAP
figure;

imagesc(-map);
set(gca, 'YTick', 0.5:(L+0.5), 'YTickLabel', floor(low_f*10)/10);
set(gca, 'XTick', J/T*(0:T) + 0.5);
set(gca, 'XTickLabel', (get(gca, 'XTick')-0.5)/J*T);
axis xy;


colormap pink;
grid on;
if PRINT_FIG
    export_fig([figpath 'video_mrdmd_map'],'-png','-transparent');
end

%% PLOT TRUE SPATIOTEMPORAL DYNAMICS

figure;

for i=1:length(modes)
    mode = -modes(i).u;
    %minmax normalization
    mode = (mode - min(mode(:)))./range(mode(:));
    figure(gcf); h=surfc(reshape(mode,nx,ny)); view(2), shading interp,
    cmap=colormap(bone(32)); brighten(-0.7);
    axis square, axis off, hold on 
    
    fname = [figpath 'video_mode' num2str(i)]; 
    if PRINT_FIG
        export_fig(fname,'-png','-r600','-transparent');
    end

    delete(h(1));
    arr = linspace(0,1,16);
    h(2).LevelList = arr(1:end-1);
    h(2).LineColor = 'k'; 
 
end

[ix,iy] = ind2sub([nx ny],qdeim);
plot(ix,iy,'r*')
axis square, axis off, box on
if PRINT_FIG
    export_fig([figpath 'video_sensors'],'-pdf');
end

%% PLOT MRDMD MODES
figure;

for i=1:3
    if (i==2)
        sgn = -1; % display mode 2 negated for consistent colormap
    else 
        sgn = 1;
    end
    mode = sgn*Mr(:,i);
    %minmax normalization for color consistency across plots
    mode = (mode - min(mode(:)))./range(mode(:));
    figure(gcf); h=surfc(reshape(mode,nx,ny)); view(2), shading interp,
    cmap=colormap(bone(32)); brighten(-0.7);
    axis square, axis off, hold on 
    
    box on
    fname = [figpath 'video_mr_mode' num2str(i)];
    if PRINT_FIG
        export_fig(fname,'-png','-r600','-transparent');
    end

    delete(h(1));
    h(2).LevelList =sort(linspace(0,1,16));    
    h(2).LineColor = 'k';
    
end

[ix,iy] = ind2sub([nx ny],qdeim);
plot(ix,iy,'r*')
axis square, axis off, box on
if PRINT_FIG
    export_fig([figpath 'video_qdmd_sensors'],'-pdf');
end
%%  PLOT POD MODES

figure;

for i=1:3
    if (i==3)
        sgn = -1; % display mode 3 negated for consistent colormap
        
    else 
        sgn = 1;
    end
    mode = sgn*Ur(:,i);
    %minmax normalization
    mode = (mode - min(mode(:)))./range(mode(:));
    figure(gcf); h=surfc(reshape(mode,nx,ny)); view(2), shading interp,
    cmap=colormap(bone(32)); brighten(-0.7);
    axis square, axis off, hold on 
    
    fname = [figpath 'video_pod_mode' num2str(i)]; 
    if PRINT_FIG
        export_fig(fname,'-png','-r600','-transparent');
    end

    delete(h(1));
    h(2).LevelList =sort(linspace(0,1,16));    
    h(2).LineColor = 'k';
    
end

[ix,iy] = ind2sub([nx ny],qdeim);
plot(ix,iy,'r*')
axis square, axis off, box on
if PRINT_FIG
    export_fig([figpath 'video_qdeim_sensors'],'-pdf');
end

%%  PLOT DMD MODES


Xaug = [Xclean(:, 1:end-1); Xclean(:, 2:end)];

X = Xaug(:, 1:end-1);
Xp = Xaug(:, 2:end);

[U,S,V] = svd(X,'econ');
r = 3;
Ur = U(:,1:r); Sr = S(1:r,1:r); Vr = V(:,1:r);

Atilde = Ur' * Xp * Vr / Sr;
[W, D] = eig(Atilde);
lambda = diag(D);

Phi = Xp * Vr / Sr * (W/D);
Phi = real(Phi(1:nx*ny,:));

figure;
for i=1:3
    if (i==3)
        sgn = 1;
    else
        sgn = -1;
    end
    mode = sgn*Phi(:,i);
    %minmax normalization
    mode = (mode - min(mode(:)))./range(mode(:));
    figure(gcf); h=surfc(reshape(mode,nx,ny)); view(2), shading interp,
    cmap=colormap(bone(32));
    axis square, axis off, hold on 
    
    fname = [figpath 'video_dmd_mode' num2str(i)]; 
    if PRINT_FIG
        export_fig(fname,'-png','-r600','-transparent');
    end

    delete(h(1));
    h(2).LevelList =sort(linspace(0,1,16));    
    h(2).LineColor = 'k';
    
end

[ix,iy] = ind2sub([nx ny],qdeim);
plot(ix,iy,'r*')
axis square, axis off, box on
if PRINT_FIG
    export_fig([figpath 'video_qdeim_sensors'],'-pdf');
end

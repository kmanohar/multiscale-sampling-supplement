function display_sensors( mask, indS, cmap )
%display_sensors Displays NOAA sea surface temperature sensors on white
%with gray continents
%
%   INPUTS
%   mask: binary mask for displaying continents
%   indS: sensor location indices
%   cmap: colormap for sensors
%
%   Modified 2018/12/31

    sensors = zeros(360,180);
    
    % number sensors by sensor order (if desired, color by order can be
    % implemented)
    % This is necessary for unique() to work properly
    P = zeros(size(mask(mask==1))); 
    P(indS)=1:length(indS);    
    
    sensors(mask==1) = P;
    b = imagesc(mask');
    shading flat, colormap(gray), drawnow
    alpha 0.3
    hold on

    
    S = reshape(sensors',360*180,1);
    [~,IC,~] = unique(S);
    
    % IC(2:end) contains linear indices of sensor locations
    [I,J] = ind2sub(size(sensors'),IC(2:end));
    
    % set sensor markersize to 150 (large)
    scatter(J,I,150,cmap,'filled','MarkerEdgeColor','k');
    %alpha 0.1
    
    axis off; hold off
end


function display_fig_sst( x, mask, pivot, cbounds )
%display_fig_sst Displays single snapshot of NOAA sea surface temperature
%(SST) data, and sensors if specified
%
%   INPUTS
%   x: single snapshot of NOAA SST temperature field
%   mask: binary mask for displaying continents
%   pivot: sensor locations (not displayed if empty)
%   cbounds: = [CLOW CHIGH] can specify the color scaling
%
%   Modified 2018/12/31

    snapshot = NaN*zeros(360*180,1);
    snapshot(mask==1) = x;
    
    sensors = zeros(360,180);
    P = zeros(size(x)); P(pivot)=ones(size(P(pivot)));
    sensors(mask==1) = P;
    
    C = reshape(real(snapshot),360,180)';
    if (~isempty(cbounds))
        b = imagesc(C,cbounds);
    else 
        b = imagesc(C);%,[-1.799 30.77999]);
    end
    
    set(b,'AlphaData',(~isnan(C)));
    if (~isempty(sensors))
        hold on
        label = linspace(min(x), max(x), length(pivot));
        S = reshape(sensors,360,180)';
        [I,J] = find(S>0);
        
        scatter3(J,I,10*ones(size(pivot)),40,label,'o','filled',...
        'MarkerEdgeColor','black');
        hold off
    end
    axis off
    shading interp, colormap(jet), drawnow
end


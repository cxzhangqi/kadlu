function [] = plot_bathy(x,y,z,CAXIS,ifig)

if nargin == 5,
    figure(ifig);
end

if ~exist('CAXIS','var')||isempty(CAXIS),
    CAXIS = [min(z(:)) max(z(:))];
end

surf(x,y,abs(z));

shading interp
ch = colorbar; set(get(ch,'ylabel'),'string','Depth (m)'); set(ch,'tickdir','out')

% caxis([min(z(:)) max(z(:))])
caxis(CAXIS)
colormap(flipud(summer));

lighting('gouraud');
camlight(40,-10)
view(2);  % These angle/elevation values seem to show 
              % the topography better for this particular dataset
hl = light;

return

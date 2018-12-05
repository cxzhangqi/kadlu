addpath('./subroutines/');
if ~exist('PE','var')||isempty(PE),
    matdir = './mat/';
    matfile = 'BaffinBay_run2_cyliNx2D__200Hz_phone274.50m.mat';%'BaffinBay_run1_RCV50m_cyli3D__200Hz_phone50.00m.mat'; 
    PE = load([matdir matfile]);
    spnglyr_ratio = 4;
    AXIS = [-65 -60 74.5 76];%[-60 -56 73 74.5]; %AXIS = [-68 -55 74 76.5];
    Bathy_dx = 50;  % bathymetry grid size
    Bathy_dy = 50;
    Bathy_Medfilt2Size = [250 250];  % use the median filter to smooth the bathymetry data (in meters)
    TriScatteredInterp_method = 'linear';
    bathy = load('../bathy/BaffinBay_bathy.mat');  
    Bathy = sub_create_BathyGridforplot(bathy.lar,bathy.lor,bathy.zr,Bathy_dx,Bathy_dy,...
        Bathy_Medfilt2Size,TriScatteredInterp_method,...
        AXIS(1:2),AXIS(3:4),PE.ENV.lon_ref,PE.ENV.lat_ref);
    [Bathy.LON,Bathy.LAT] = meshgrid(Bathy.x,Bathy.y);
    [Bathy.LON(:),Bathy.LAT(:)] = sub_transfer_XY_to_LL(Bathy.LON(:),Bathy.LAT(:),PE.ENV.lon_ref,PE.ENV.lat_ref);
end

RCV_depth = 50; 

    % 3D TL solution
    [~,idz] = min(abs(PE.Ez_z-RCV_depth));
    RCV_depth = PE.Ez_z(idz); 
    tl0 = fftshift(-20*log10(squeeze(abs(PE.Ez_nx2d(idz,:,:)))),1);
    theta0 = fftshift(PE.Ez_theta_nx2d,1); ind_half = size(PE.Ez_theta_nx2d,1)/2;
    y0 = min(theta0(:,end));%/spnglyr_ratio*(spnglyr_ratio-1);
    y1 = max(theta0(:,end));%/spnglyr_ratio*(spnglyr_ratio-1);
    
    ind = theta0>=y0 & theta0<=y1; 
    
    THETA = theta0(ind);
    R = PE.r_nx2d;
    TL = tl0(ind,:); 
    
    [R,THETA] = meshgrid(R,THETA);
    [X,Y] = pol2cart(THETA,R); LON = nan*X; LAT = nan*Y;
    [LON(:), LAT(:)] = sub_transfer_XY_to_LL(X(:),Y(:),PE.ENV.lon_ref,PE.ENV.lat_ref,PE.ENV.rot_ref);

figure(12)
% pcolor(LON,LAT,TL-10*log10(ones(size(TL,1),1)*PE.r)); hold on; shading flat;
pcolor(LON,LAT,TL); hold on; shading flat;
[c,h]=contour(Bathy.LON,Bathy.LAT,Bathy.WD,[200 300 500 750],'linecolor',[1 1 1]*0);
% caxis([15 35]+3);
set(gca,'DataaspectRatio',[1 cos(mean(LAT(:))/180*pi) 1],'tickdir','in')
axis(AXIS); %([-68 -55 74 76.5])
% set(gca,'xtick',[122.5:.025:122.65],'ytick',[25.5:.025:25.65])
ih=colorbar('ydir','reverse'); %ih=colorbar('ydir','reverse','position',[.334 .275 .02 .25],'tickdir','out','ytick',[18:4:38]); 
set(get(ih,'title'),'string','dB')
h1 = clabel(c,h,'labelspacing',150);
for itmp = 1:length(h1); set(h1(itmp),'string',[get(h1(itmp),'string') ' m'],'color',[1 1 1]*.5,'fontsize',8); end
% title({sprintf('TL at %.2f m depth',RCV_depth) '(cylindrical spreading loss removed)'},'fontsize',12)
title(sprintf('TL from sources at %.2f m depth to a receiver at 274 m',RCV_depth),'fontsize',12)
% title(sprintf('TL at %.2f m depth',RCV_depth),'fontsize',12)
xlabel('Longitude','fontsize',12); ylabel('Latitude','fontsize',12)


    
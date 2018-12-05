clear
global WD SSP
addpath subroutines
addpath /Users/dbarclay/Documents/MatlabScripts/Bathymetry/map

lon_ref = 142.1;
lat_ref = 11.2;
rot_ref = 0;

% bathymetry
bathy = load('env_database/BathyData_Mariana_500kmx500km.mat');
bathy.longrat = bathy.longrat-360;
[X,Y] = sub_transfer_LL_to_XY(bathy.longrat,bathy.latgrat,lon_ref,lat_ref,rot_ref);

F = scatteredInterpolant(X(:), Y(:), -bathy.mat(:),'natural');

xkm=deg2km(bathy.longrat);
ykm=deg2km(bathy.latgrat);
xkm=xkm-mean(xkm(1,:));
ykm=ykm-mean(ykm(:,1));



% % for rot_ref = 0; 
if(rot_ref==0) 
    dx = fix(mean(mean(diff(X,1,2)))/100)*100;
    dy = fix(mean(mean(diff(Y,1,1)))/100)*100; dy = -dy;
    x = dx*(fix(max(X(:,1))/dx):fix(min(X(:,end))/dx));
    y = dy*(fix(max(Y(end,:))/dy):fix(min(Y(1,:))/dy)).';
elseif(rot_ref==90)
    dx = fix(mean(mean(diff(X,1,1)))/100)*100; dx = -dx;
    dy = fix(mean(mean(diff(Y,1,2)))/100)*100; dy = -dy; 
    x = dx*(fix(max(X(end,:))/dx):fix(min(X(1,:))/dx));
    y = dy*(fix(max(Y(:,end))/dy):fix(min(Y(:,1))/dy)).';
end

[X,Y] = meshgrid(x,y);

WD.field = F(X,Y).';


WD.dx = dx; WD.x0 = x(1); WD.x1 = x(end);
WD.dy = dy; WD.y0 = y(1); WD.y1 = y(end);
WD.lon_ref = lon_ref; WD.lat_ref = lat_ref; WD.rot_ref = rot_ref;

figure;
subplot(221); imagesc(bathy.longrat(1,:),bathy.latgrat(:,1),-bathy.mat); colorbar; axis xy
subplot(223); imagesc(x,y,WD.field.'); colorbar; axis xy

% water sound speed
water_profiles = load('env_database/WaterColumnProfiles_Mariana_Fake.mat');
water_profiles.Depth = water_profiles.Depth(1:34077);
water_profiles.c = water_profiles.c(1:34077);
dz = 1; dz_smooth = 10; 
z0 = 0; 
z1 = fix(max(WD.field(:))/dz)*dz; 
z = (fix(z0/dz):fix(z1/dz))*dz; z = z(:); 
ssp = nan(size(z)); 
for itmp = 1:length(z)
    if z(itmp)>=min(water_profiles.Depth(:))+dz_smooth/2 &&  z(itmp)<=max(water_profiles.Depth(:))-dz_smooth/2,
        idz = water_profiles.Depth>z(itmp)-dz_smooth/2 & water_profiles.Depth<=z(itmp)+dz_smooth/2;
        ssp(itmp) = mean(water_profiles.c(idz));
    end
end
ssp = interp1(z(~isnan(ssp)),ssp(~isnan(ssp)),z,'linear','extrap'); 

subplot(122); plot(ssp,z,'.-',water_profiles.c,water_profiles.Depth,'r'); axis ij

SSP.field = shiftdim(ssp,-2);
SSP.dx = 0;  SSP.x0 = 0;    SSP.x1 = 0;
SSP.dy = 0;  SSP.y0 = 0;    SSP.y1 = 0;
SSP.dz = dz; SSP.z0 = z(1); SSP.z1 = z(end);

save Mariana_ENV WD SSP

rmpath subroutines

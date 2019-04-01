clear
clear global
addpath subroutines

global Ez Ez_z Ez_y  % these variabiles defines quick outputs of 3D field
Ez_z = [.1];

% --- env configuration  ---
global ENV
ENV.c0 = 1500;  % reference sound speed

load('Mariana_ENV.mat','SSP')
mm=(1:10673)/10673+1;
mm=flipdim(mm,2);
SSP.field(1,1,:)=SSP.field(1,1,:).*shiftdim(mm,-1);

NSQ.field = (ENV.c0./SSP.field).^2;
NSQ.dx = SSP.dx; NSQ.x0 = SSP.x0; NSQ.x1 = SSP.x1;
NSQ.dy = SSP.dy; NSQ.y0 = SSP.y0; NSQ.y1 = SSP.y1;
NSQ.dz = SSP.dz; NSQ.z0 = SSP.z0; NSQ.z1 = SSP.z1;
ENV.NSQ = NSQ; clear NSQ
ENV.ndx_ChangeNSQ = inf;  % range-independent ssp

load('Mariana_ENV.mat','WD')
ENV.WD = WD; clear WD
ENV.ndx_ChangeWD = 3; %1;  % how often to update the env

ENV.rhow = 1;         % water density
ENV.cb = 1700;        % homogeneous bottom sound speed
ENV.bloss = 0.5 ;     % homogeneous bottom attenuation db/lambda
ENV.rhob = 1.5;       % homogeneous bottom density

smoothing_length_rho = []; % smoothing length on the water/bottom interface
smoothing_length_ssp = []; % empty uses default

% -------------------------------------
icase = 'Mariana_Nx2D_longrng';
isSingle = 1; isplot = 0;

% --- acoustic source  ---
xs = 0;
ys = 0;
ZS = [sub_SeafloorDepth(xs,ys)-95; sub_SeafloorDepth(xs,ys)-100; sub_SeafloorDepth(xs,ys)-105];   % m
izs = input(sprintf('which zs? (1) %.2f m (2) %.2f m (3) %.2fm  ',ZS(1),ZS(2),ZS(3))); 
zs = ZS(izs); 
freq = [10];  % Hz

% --- PE calculation parameters  ---
c0 = ENV.c0;
lambda0 = c0./freq(end);  % reference wavelength, choose smallest one

% r  (PE marching direction)
%%%OLI steplength = lambda0/2;  % dx in meters
steplength = 1000;  % dx in meters
rmax = 50e3;% m
ndxout = 1; %2; 
numstep = round(rmax/steplength);  % number of marching steps and distance

% theta
Ltheta = 2*pi;%*model_domain.ThinknessOfArtificialAbsorpLayer_ratio_y/(-1+model_domain.ThinknessOfArtificialAbsorpLayer_ratio_y);
%%%OLI dtheta = 1/180*pi;
dtheta = 45/180*pi;
ntheta = ceil(Ltheta/dtheta); if mod(ntheta,2)==1, ntheta = ntheta+1; end
ndy_3DSliceout=1;

% z
%%%OLI dz = 10;
dz = 1000;
ThinknessOfArtificialAbsorpLayer_ratio_z = 6;   % The default value (4) will be used if it is an empty variable
nz = 2*12e3/ThinknessOfArtificialAbsorpLayer_ratio_z*(ThinknessOfArtificialAbsorpLayer_ratio_z+1)/dz;
nz = round(nz/2)*2;  % nz must be an even number
ndz_3DSliceout = 1;

% PE starter type
PEStarter.type = 'Thomson''s';
PEStarter.aperature = 88; %deg

% output settings
matdir = sprintf('mat_%s_%dHz_zs%dm',icase,freq(1),fix(zs));  % create a directory for saving mat files
if ~exist(matdir,'dir'),
    mkdir(matdir);
end

YZSlice_output_folder = []; % sprintf('%s/YZSlice',matdir);  % create a directory for saving mat files
if ~isempty(YZSlice_output_folder)&&~exist(YZSlice_output_folder,'dir'),
    mkdir(YZSlice_output_folder);   % if YZSlice_output_folder is an empty variabile, no output slide will be saved.
end

sub_EnvInput = 'default';  % env loaders
sub_Output = 'default';
if strcmpi(sub_EnvInput,'default'), clear sub_envInput_Nx2D;       % reset the functions
else eval(sprintf('clear %s;',sub_EnvInput)); end
if strcmpi(sub_Output,'default'), clear sub_output_Nx2D;
else eval(sprintf('clear %s;',sub_Output)); end

% -----------------------
% Run PE !!
% -----------------------
outfile  = sprintf('%s_%dHz.mat',icase,freq(1));

if 1%~exist([matdir '/' outfile],'file'),
    [r,theta,z,psifinal,Af] = propNx2DWAPE(...
        icase,...
        freq,xs,ys,zs,...
        PEStarter,c0,...
        numstep,ntheta,nz,...
        steplength,dtheta,dz,...
        sub_EnvInput,sub_Output,...
        smoothing_length_rho,smoothing_length_ssp,...
        ThinknessOfArtificialAbsorpLayer_ratio_z,...
        isSingle,...
        ndxout,ndy_3DSliceout,ndz_3DSliceout,...
        YZSlice_output_folder,isplot);
    save([matdir '/' outfile],'r','theta','z','psifinal','Af','Ez','Ez_z','Ez_y','freq','ENV') %'-v7.3')
else
    load([matdir '/' outfile],'r','theta','z','psifinal','Af','Ez','Ez_z','Ez_y','freq','ENV')
end

plot_r = r(2:end);
plot_theta = fftshift(squeeze(Ez_y(:,:,1)));
[R,TH] = meshgrid(plot_r,plot_theta);
[X,Y] = pol2cart(TH,R);


for idz = 1:length(Ez_z);
    
    plotfreq=1;
    
    ifig = figure; clf
    set(gcf,'papersize',[8.5 11],'paperposition',[.25 .25 8 10.5])
    
    SPfield = fftshift(squeeze(Ez(idz,:,2:end,plotfreq)),1);
    pcolor(X/1e3,Y/1e3,(20*log10(abs(SPfield)))); shading flat
    
    caxis([-60 0]+max(caxis))
    ih = colorbar; set(get(ih,'ylabel'),'string','dB')
    xlabel('X (km)');ylabel('Y (km)')
    title(sprintf('Nx2D PE solution at Z = %.2f m, frequency %d Hz',Ez_z(idz),freq))
    
    % set(gca,'position',[.13 .11 .65 .815])
    set(gca,'tickdir','out'); % axis equal; axis tight; axis([-1 1 -1 1]*max(plot_r(1,:))/1e3);
    axis equal; axis xy; axis tight
    
    picfile = sprintf('Nx2DPEsolutionatZ%dm.png',fix(Ez_z(idz)));
    print('-dpng','-r300',[matdir '/' picfile])
end

figure;
imagesc((ENV.WD.x0:ENV.WD.dx:ENV.WD.x1)/1e3,(ENV.WD.y0:ENV.WD.dy:ENV.WD.y1)/1e3,ENV.WD.field); colorbar; axis equal; xlabel('x (km)'); ylabel('y (km)'); title('Bathymetry')
hold on; plot((rmax*cos((0:359)/180*pi)+xs)/1e3,(rmax*sin((0:359)/180*pi)+ys)/1e3,'w-','linewidth',2); axis xy
picfile = sprintf('bathymetry.png');
print('-dpng','-r300',[matdir '/' picfile])

figure;
imagesc(r/1e3,z,20*log10(abs(Af(:,:,plotfreq)))); colorbar; axis ij
xlabel('r (km)'); zlabel('z (m)')
hold on; plot(r/1e3,sub_SeafloorDepth(r(:),zeros(size(r(:)))),'w','linewidth',2)
picfile = sprintf('Nx2DPEsolutionalongtheta0deg.png');
print('-dpng','-r300',[matdir '/' picfile])

rmpath subroutines


function [x,y,z,psifinal,Af] = propNx2DWAPE(...
    icase,...
    freq,xs,ys,zs,...
    PEStarter,c0,...
    numstep,ny,nz,...
    steplength,dy,dz,...
    sub_EnvInput,sub_Output,...
    smoothing_length_rho,smoothing_length_ssp,...
    ThinknessOfArtificialAbsorpLayer_ratio_z,...
    isSingle,...
    ndxout,ndy_3DSliceout,ndz_3DSliceout,...
    YZSlice_output_folder,isplot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Require even number of grid points
% Model grid is defined in a way like y = [0 1 2 -3 -2 -1]
% y is the azimuth variable
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% declare PERSISTENT variables
global Y Z k0 kz nfft
global n2in denin ddenin d2denin
global Ez Ez_z Ez_y wd_mask

if isSingle,
    disp('----- Nx2D Standard Wide Angle SSF PE (single precision) -----')
else
    disp('----- Nx2D Standard Wide Angle SSF PE (double precision) -----')
end

if ~exist('ThinknessOfArtificialAbsorpLayer_ratio_z','var')||isempty(ThinknessOfArtificialAbsorpLayer_ratio_z),
    ThinknessOfArtificialAbsorpLayer_ratio_z = 4; disp('      Set ThinknessOfArtificialAbsorpLayer_ratio_z = 4'); end
if ~exist('isplot','var')||isempty(isplot), isplot = 0; end
if strcmpi(sub_EnvInput,'default'), envInput = @sub_envInput_Nx2D; 
else
    if strfind(sub_EnvInput,'.m'), sub_EnvInput = sub_EnvInput(1:strfind(sub_EnvInput,'.m')-1); end
    eval(sprintf('envInput = @%s;',sub_EnvInput))
end
if strcmpi(sub_Output,'default'), modelOutput = @sub_output_Nx2D; 
else
    if strfind(sub_Output,'.m'), sub_Output = sub_Output(1:strfind(sub_Output,'.m')-1); end
    eval(sprintf('modelOutput = @%s;',sub_Output))
end

% --- reference wavenumber ---
lambda0=c0/freq ;
k0=2*pi*freq/c0;
k0sq=k0.^2;

% --- set up grid ---
dx=steplength;
x=[0 (ndxout:ndxout:numstep)*steplength];

Ly=dy*ny;   % rough estimate of angular aperature
y=[(1:ny/2)-1 (-ny/2:-1)]*dy;

Lz=nz*dz;
z=[(1:nz/2)-1 (-nz/2:-1)]*dz; z = z(:);
kz=[(1:nz/2)-1 (-nz/2:-1)]*2*pi/Lz; kz = kz(:);

% --- initialize the grid ---
Ez_y = nan(ny,1);
Ez = nan(length(Ez_z),ny,length(x));
[Y,Z]=meshgrid(y,z);
if isSingle,
    psi = zeros(size(Z),'single');    % allocate memory for the following variables
    fr_half = zeros(size(Z),'single');
    fr_full = zeros(size(Z),'single');
    atten0 = nan(size(Z),'single');
    U = zeros(size(Z),'single');
    Af=nan(nz/2,length(x),'single');
else
    psi = zeros(size(Z));    % allocate memory for the following variables
    fr_half = zeros(size(Z));
    fr_full = zeros(size(Z));
    fr_Nx2D = zeros(nz,1);
    atten0 = nan(size(Z));
    U = zeros(size(Z));
    Af=nan(nz/2,length(x));
end
nfft = 0;

% --- boundary conditions on z --- Note that sponge layer on y varys when the angular interval angles
% using sponge layers to simulate radiation boundary condition
% see Section 6.5.3 of Computational Ocean Acoustics for the formula using here,
% i.e.,    n^2 = nb^2 + i * alpha * exp( -(z-zmax).^2/D.^2 )
ArtificialAbsorpCoeff =  1/log10(exp(1))./lambda0*2/k0;     % 20dB to zmax and ymax, Gaussian decayed length: D

ThinknessOfArtificialAbsorpLayer = max(abs(z))/ThinknessOfArtificialAbsorpLayer_ratio_z;      % which gives effective PE domain within zmax-ThinknessOfArtificialAbsorpLayer
D = ThinknessOfArtificialAbsorpLayer/3;         % Gaussian decayed length: D
atten0(:) = sqrt(-1) * ArtificialAbsorpCoeff * exp(-(abs(Z(:))-max(abs(z))).^2/D.^2);      % the sponge layer on z
clear ArtificialAbsorpCoeff ThinknessOfArtificialAbsorpLayer D

if isplot, % plot the sponge layer
    figure(37); atten_temp=imag(k0*sqrt(1+atten0))*20*(c0/freq)*log10(exp(1));
    imagesc(fftshift(y),fftshift(z),fftshift(atten_temp)); clear atten_temp;
    xlabel('Y');ylabel('Z');title('Sponge boundary');
    ih=colorbar; set(get(ih,'ylabel'),'string','dB/\lambda');drawnow
end

% --- initial condition, i.e., PE Starter ---
psi(:,:) = sub_CylindPEStarter(k0,kz,PEStarter,zs)*ones(1,ny);

% --- initial Nx2D free propagator
fr_half(:,:)=exp(1i*dx/2*( sqrt(k0sq-kz.^2) - k0 ))*ones(1,ny);
fr_full(:,:)=exp(1i*dx*( sqrt(k0sq-kz.^2) - k0 ))*ones(1,ny);

% --- Default values of smoothing_length_rho and smoothing_length_ssp ---
if ~exist('smoothing_length_rho','var')||isempty(smoothing_length_rho),
    smoothing_length_rho = c0/freq/4;
    disp('      smoothing_length_rho = lambda0/4')
end
if ~exist('smoothing_length_ssp','var')||isempty(smoothing_length_ssp),
    smoothing_length_ssp = eps;
    disp('      smoothing_length_ssp = eps')
end

% output field at 0
dista = 0;
iNextOutput = 1;
Af(:,iNextOutput) = modelOutput(icase,psi,dista,ndy_3DSliceout,ndz_3DSliceout,YZSlice_output_folder,isplot);
iNextOutput = iNextOutput+1;

% --- PE marching starts here
is_halfstep = 1;    % initially, half step to dx/2
tic; seco=toc; 
for jj=1:numstep
    
    %(1) x --> x+dx/2 free propagation
    if is_halfstep, psi = fr_half.*psi; end
    
    %(2) do phase adjustment at x+dx/2
    [isnewscreen,isupdate] = envInput(freq,c0,dista+dx/2,dx,xs,ys,smoothing_length_rho,smoothing_length_ssp);
    if isnewscreen,
        %   generate phase screen for this range
        %   generate phase screen for this range
        % U = sqrt( n2in + atten0 + 1/2/k0/k0*( (rhob-rhow)*d2Hin./denin - 3/2*((rhob-rhow)*dHin./denin).^2 ) ) - 1;
        % %         ^^^^   ^^^^^^   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        % %         n^2  sponge layer                            density effect
        % U = exp(1i*dx*k0*U);    % phase screen
        U(:,isupdate) = exp(1i*dx*k0*(-1 + sqrt( n2in(:,isupdate) + atten0(:,isupdate) + ...
            1/2/k0/k0*(d2denin(:,isupdate)./denin(:,isupdate) - 3/2*(ddenin(:,isupdate)./denin(:,isupdate)).^2) ) ) );
    end
    psi = fft(U.*ifft(psi));    nfft = nfft+2; 
    
    %(3) x+dx/2 --> x+dx free propagation
    dista = dista + dx;
    % output field?
    if any(x==dista),
        is_halfstep = 1;        % output field at x+dx, so half step from x+dx/2
        psi = fr_half.*psi;
        Af(:,iNextOutput) = modelOutput(icase,psi,dista,ndy_3DSliceout,ndz_3DSliceout,YZSlice_output_folder,isplot);
        iNextOutput = iNextOutput+1;
    else                        % if not output filed, full step to x+dx/2+dx
        is_halfstep = 0;
        psi = fr_full.*psi;
    end
    
    if (toc-seco)>10, 
        seco=toc; 
        fprintf('    Elapsed time to %d/%d is %.2f sec to %.2f m with %d ffts.\n',jj,numstep,seco,dista,nfft)
    end
    
end

psifinal = ifft(psi)*exp(1i*k0*dista)/sqrt(dista).*sqrt(denin);   nfft = nfft+1; 
psifinal = fftshift(psifinal(1:nz/2,:),2);
z=z(1:nz/2);
y=fftshift(y);
fprintf('----- Finsih in %.2f sec for a %.2f m run with %d ffts. -----\n\n',toc,dista,nfft)
    
clear global Y Z k0 kz nfft
clear global n2in denin ddenin d2denin
clear global wd_mask

if strcmpi(sub_EnvInput,'default'), clear sub_envInput_Nx2D;
    else eval(sprintf('clear %s;',sub_EnvInput)); end
if strcmpi(sub_Output,'default'), clear sub_output_Nx2D;
    else eval(sprintf('clear %s;',sub_Output)); end

return


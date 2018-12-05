function [isnewnsq,isnew_env] = sub_envInput_Nx2D(freq,c0,dista,dx,xs,ys,smoothing_length_rho,smoothing_length_ssp)

% use the GLOBAL variable declared before
global ENV
global Y Z
global n2in denin ddenin d2denin 
global wd_mask

% declare PERSISTENT variables seen only in this sub_envInput subroutines
persistent wd_x_next 
persistent wd_new wd_old Z_sub_wd 
persistent n2b rhob H_c H_rho DwdDy_new DwdDy % DDwdDyy_new DDwdDyy
persistent NSQ_x_next 
persistent n2w_new n2w rhow 
persistent costheta sintheta

nf=length(freq);

if isempty(wd_new),   % allocate variables for faster env updates
    wd_new = nan(1,size(Y,2));   wd_old = nan(1,size(Y,2));
    DwdDy_new = nan(1,size(Y,2));   DwdDy = nan(size(Z));
    % DDwdDyy_new = nan(1,size(Y,2));   DDwdDyy = nan(size(Z));
    Z_sub_wd = nan([size(Z) nf]);   wd_mask = nan([size(Z) nf]);
    H_rho = nan([size(Z) nf]);   ddenin = nan([size(Z) nf]);   d2denin = nan([size(Z) nf]);
    denin = nan([size(Z) nf]);
    wd_x_next = 0;  NSQ_x_next = 0;
    H_c = nan([size(Z) nf]);
    n2w_new = zeros([size(Z) nf]);  n2w = zeros([size(Z) nf]);
    n2in = nan([size(Z) nf]);
    costheta = cos(Y(1,:)); sintheta = sin(Y(1,:)); 
end

% initiate bottom model
if (dista==dx/2),
    cb = ENV.cb;        % bottom sound speed
    bloss = ENV.bloss;    % bottom attenuation db/lambda
    rhob = ENV.rhob;       % bottom density
    
    % complex bottom sound speed
    kbi = bloss./(cb./freq)/20/log10(exp(1));
    betab = kbi/2/pi./freq;
    for(kk=1:length(freq))
        cbi(kk,:) = roots([betab(kk),-1,betab(kk)*cb.^2]);
    end
    cbi = cbi(imag(cbi)==0); cbi = cbi(cbi>=0&cbi<cb);
    cb = cb - 1i*cbi;
    n2b = (c0./cb).^2;
end

x = xs+costheta*dista; y = ys+sintheta*dista;

isnewnsq = 0;
% bathymetry
isnewbathy = false(1,size(Y,2));
if (dista==dx/2)||(dista>=wd_x_next),
    wd_x_next = dista + ENV.ndx_ChangeWD*dx;
    % get bathymetry
    [wd_new(:),DwdDy_new(:)]  = sub_SeafloorDepth(x,y);
    isnewbathy = (wd_old~=wd_new)&(~isnan(wd_new));
    wd_old = wd_new;
    wd_mask(:,isnewbathy) = ones(size(Z,1),1)*wd_new(isnewbathy); % water depth mask
    DwdDy(:,isnewbathy) = ones(size(Z,1),1)*DwdDy_new(isnewbathy); %
    % DDwdDyy(:,isnewbathy) = ones(size(Z,1),1)*DDwdDyy_new(isnewbathy); %
    Z_sub_wd(:,isnewbathy) = abs(Z(:,isnewbathy))-wd_mask(:,isnewbathy);
    % fprintf('   Updating bathymetry at %.2f m\n',dista)
end

% water column
isnewNSQ=false(1,size(Y,2));

[nx,ny,nothing]=size(ENV.NSQ.field);
if nx == 1 && ny == 1, range_independent = 1; else range_independent = 0; end

if (dista==dx/2)||((dista>=NSQ_x_next)&&~range_independent),  % update water column
    NSQ_x_next = dista + ENV.ndx_ChangeNSQ*dx;
    % interpolate water nsq
    n2w_new = n2w;
    IDZ = (Z<=0) ;%yt01152013& (Z_sub_wd<=200);
    [nothing,idy] = find(IDZ); IDZ = find(IDZ(:));
    n2w_new(IDZ) = sub_NSQ(x(idy).',y(idy).',Z(IDZ));
    clear IDZ
    isnewNSQ = any((n2w([1 size(Z,1):-1:size(Z,1)/2+2],:) ...
        -n2w_new([1 size(Z,1):-1:size(Z,1)/2+2],:))~=0,1);
    n2w_new(2:size(Z,1)/2,isnewNSQ) = n2w_new(size(Z,1):-1:size(Z,1)/2+2,isnewNSQ);
    n2w(:,isnewNSQ) = n2w_new(:,isnewNSQ);
    rhow = ENV.rhow;       % water density
    % fprintf('   Updating water column at %.2f m\n',dista)
end

isnew_env = isnewbathy | isnewNSQ;

if ~isempty(find(isnew_env,1)),
    % smooth ssp
    H_c(:,isnew_env) = (1 + tanh(Z_sub_wd(:,isnew_env)/smoothing_length_ssp/2) )/2;
    n2in(:,isnew_env) = n2w(:,isnew_env)+(n2b-n2w(:,isnew_env)).*H_c(:,isnew_env);
    itmp = (wd_mask(1,:)==0); if any(itmp), n2in(1,itmp) = n2b; end
    
    % smooth density
    TANH = tanh(Z_sub_wd(:,isnew_env)/smoothing_length_rho/2);
    H_rho(:,isnew_env) = (1 + TANH) /2;
    denin(:,isnew_env) = rhow+(rhob-rhow)*H_rho(:,isnew_env);
    
    SECH2 = 1./cosh(Z_sub_wd(:,isnew_env)/smoothing_length_rho/2);
    SECH2 = SECH2 .* SECH2;
    ddenin(:,isnew_env) =  SECH2 /smoothing_length_rho/2 .* sqrt(1+DwdDy(:,isnew_env).^2);
    ddenin(:,isnew_env) =  (rhob-rhow)/2 * ddenin(:,isnew_env);
    %                                                                  ----------------------------------------------------
    d2denin(:,isnew_env) = - SECH2 /smoothing_length_rho/2 .* ( TANH/smoothing_length_rho.*(1+DwdDy(:,isnew_env).^2) );%+ DDwdDyy(:,isnew_env) );
    d2denin(:,isnew_env) =  (rhob-rhow)/2 * d2denin(:,isnew_env);
    
    % fprintf('   Updating phase screen at %.2f m\n\n',dista)
    isnewnsq = 1;
end

return


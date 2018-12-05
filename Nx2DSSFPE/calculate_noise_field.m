addpath subroutines

global ENV
load('Tonga_ENV.mat','WD')
ENV.WD = WD; clear WD

icase = 'Tonga_Nx2D_longrng';   % name of the study case

xs = 0;
ys = 0;
zs = [sub_SeafloorDepth(xs,ys)-95; sub_SeafloorDepth(xs,ys)-105; sub_SeafloorDepth(xs,ys)-100];   % m

freq = 50;  % Hz

% load PE outputs for different PE starter depths
if ~exist('Output','var'),
    for izs = 1:length(zs);
        matdir = sprintf('mat_%s_%dHz_zs%dm',icase,freq,fix(zs(izs)));  % create a directory for saving mat files
        outfile  = sprintf('%s_%dHz.mat',icase,freq);
        Output(izs) = load([matdir '/' outfile],'r','theta','z','psifinal','Af','Ez','Ez_z','Ez_y');
    end
end
 
% compute the normalization factor (dr and dtheta) 
for izs = 1:length(Output);
    Output(izs).zs = zs(izs);
    Output(izs).dth = nanmean(diff(fftshift(Output(izs).Ez_y)));
    Output(izs).dr = nanmean(diff(Output(izs).r));
end


% compute cross spectral density matrix
dtheta = 5; theta = 0:dtheta:360; theta = theta(1:end-1).';
dtheta = dtheta/180*pi;  theta = theta/180*pi;
for izs = 1:length(Output);
    for jzs = 1:izs;
        CSD(izs,jzs).sigsq = squeeze(Output(izs).Ez(1,:,:)) .* conj(squeeze(Output(jzs).Ez(1,:,:)));
        CSD(izs,jzs).sigsq = CSD(izs,jzs).sigsq .* (Output(izs).dr*ones(size(Output(izs).Ez_y,1),1)*(Output(izs).r.*Output(izs).dth));
        CSD(izs,jzs).SIGSQ = nan(length(theta),length(Output(izs).r));
    end
end
for idr = 2:length(Output(1).r);
    Ez_theta = Output(1).Ez_y(:);
    Ez_theta(Ez_theta<0) = 2*pi+Ez_theta(Ez_theta<0);
    for ith = 1:length(theta)
        theta0 = theta(ith)-dtheta/2; theta1 = theta(ith)+dtheta/2;
        if theta0 < 0,
            theta0 = 2*pi+theta0;
            itmp = Ez_theta>theta0 | Ez_theta<=theta1;
        else
            itmp = Ez_theta>theta0 & Ez_theta<=theta1;
        end
        for izs = 1:length(Output);
            for jzs = 1:izs;
                CSD(izs,jzs).SIGSQ(ith,idr) = mean(CSD(izs,jzs).sigsq(itmp,idr))/Output(izs).dth*dtheta;
            end
        end
    end
end


for izs = 1:length(Output);
    for jzs = 1:izs;
        CSD(izs,jzs).SIGSQ(isnan(CSD(izs,jzs).SIGSQ)) = 0; 
        CSD(izs,jzs).SIGSQ_cum = cumsum(CSD(izs,jzs).SIGSQ,2); 
    end
end


for izs = 1:length(Output);
    for jzs = 1:izs;
        if izs == jzs;
            CSD(izs,jzs).coh = CSD(izs,jzs).SIGSQ_cum;
        else
            CSD(izs,jzs).coh = CSD(izs,jzs).SIGSQ_cum ./ sqrt(CSD(izs,izs).SIGSQ_cum .* CSD(jzs,jzs).SIGSQ_cum) ;
        end
    end
end



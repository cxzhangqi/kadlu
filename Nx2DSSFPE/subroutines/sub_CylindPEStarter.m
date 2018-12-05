function [psi] = sub_CylindPEStarter(k0,kz,Starter,zs)

% use global variables declared before
global Z

% the following three starters simulates point source
switch lower(Starter.type(1:2));
    case 'ga'
        %(1) Gaussian starter
        psi = sqrt(k0) * exp( -0.5*k0.^2 *(Z(:,1)-zs).^2 ); 
        psi = psi - ( sqrt(k0) * exp( -0.5*k0.^2 *(Z(:,1)+zs).^2 ));
        psi = fft(psi);
    case 'gr'
        %(2) Greene's starter 
        psi = sqrt(k0) * (1.4467-.04201*k0*k0*(Z(:,1)-zs).^2) .*exp(-(k0*k0*(Z(:,1)-zs).^2)/3.0512 ) ;
        psi = psi - ( sqrt(k0) * (1.4467-.04201*k0*k0*(Z(:,1)+zs).^2) .*exp(-(k0*k0*(Z(:,1)+zs).^2)/3.0512 ) ); %   d1 and d2 are distances already squared....
        psi = fft(psi);
    case 'th'
        %(3) Thomson's starter
        psi = exp(-1i*pi/4)*2*sqrt(2*pi)*sin(kz*zs)./sqrt(sqrt(k0.^2-kz.^2));
        % normalize the starter
        psi = psi/(Z(2,1)-Z(1,1));
        psi(size(Z,1)/2+1,:) = 0; 
        % taper the spectrum to obtain desired angle using Turkey window
        kcut1 = k0*sin(Starter.aperature/180*pi); kcut0 = k0*sin((Starter.aperature-1.5)/180*pi);
        W = .5*(1+cos(pi/(kcut1-kcut0)*(abs(kz)-kcut0)));
        W(abs(kz)>=kcut1) = 0; W(abs(kz)<=kcut0) = 1;
        psi = psi.*W; psi(abs(kz)>=kcut1) = 0;
end

return
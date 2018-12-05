function [Af] = sub_output_Nx2D(icase,psi,dista,ndy_3DSliceout,ndz_3DSliceout,YZSlice_output_folder,isplot)

% use PERSISTENT variables declared before
global Y Z k0 kz nfft
global wd_mask n2in denin
global Ez Ez_z Ez_y

% declare PERSISTENT variables seen only in this sub_output subroutines
persistent idz_out idy_out z_3DSlice y_3DSlice Ez_ifft_kernel iout_Ez

% --- which points to be output?
yso=find(Y(1,:)>=0,1);
nzhalf=size(Z,1)/2;
force_3DSliceout_x = []; 

% output 3D volumn slice
if ~isempty(YZSlice_output_folder)&&isempty(idz_out),
    idz_out = 1:ndz_3DSliceout:size(Z,1)/2;
    idy_out = 1:ndy_3DSliceout:size(Y,2)/2;
    idy_out = [fliplr(size(Y,2)+1-ndy_3DSliceout:-ndy_3DSliceout:size(Y,2)/2+1) idy_out];
    z_3DSlice = Z(idz_out,1);
    y_3DSlice = Y(1,idy_out);
end

% if isempty(idz_Ez)
%     idz_Ez = nan(size(Ez_z));
%     IDZ = find(Ez_z~=Inf);
%     for itmp = IDZ(:).';
%         [junk,idz_Ez(itmp)] = min(abs(Z(:,1)-Ez_z(itmp)));
%     end
%     Ez_z(IDZ) = Z(idz_Ez(IDZ),1);
%     idz_Ez(Ez_z==Inf) = Inf;
%     iout_Ez = 0;
%     Ez_y = Y(1,:).';
% end
if isempty(Ez_ifft_kernel)
    Ez_ifft_kernel = exp(1i*Ez_z(:)*kz(:).')/length(kz);
    iout_Ez = 0;
    Ez_y = Y(1,:).';
end

iout_Ez = iout_Ez + 1;
if dista ~= 0,
    % psifinal = ifft(psi)*exp(1i*k0*dista)/sqrt(dista).*sqrt(denin);
    % Ez((Ez_z~=Inf),:,iout_Ez) = abs(psifinal(idz_Ez(Ez_z~=Inf),:));
    % if any(Ez_z == Inf),
    %     tmp=psifinal; %depth-integrated sound intensity over the water column
    %     tmp(Z<-wd_mask&Z>0)=0;
    %     Ez((Ez_z==Inf),:,iout_Ez)=sum(abs(tmp).^2).*(Z(2)-Z(1)); %depth-integrated sound intensity over the water column
    % end
    Ez(:,:,iout_Ez) = (Ez_ifft_kernel*psi)*exp(1i*k0*dista)/sqrt(dista).*sqrt(denin(round(Ez_z/(Z(2)-Z(1)))+1,:)); 
    psi = ifft(psi)*exp(1i*k0*dista)/sqrt(dista).*sqrt(denin);   nfft = nfft+1; 
else
    psi = ifft(psi);   nfft = nfft+1; 
end
Af=psi(1:nzhalf,yso); %sound intensity on the vertical x-z plane crossing the source

% output 3D volumn slice
if ~isempty(YZSlice_output_folder),
    x_3DSlice = dista;
    outfile = sprintf('%s_YZSlices_X(%09.2fm).mat',icase,x_3DSlice);
    YZ_Slice = psi(idz_out,idy_out);
    save([YZSlice_output_folder '/' outfile],'x_3DSlice','y_3DSlice','z_3DSlice','YZ_Slice')
end

if isplot&&(mod(iout_Ez,50)==1),
    [theta,isort] = sort(Y(1,:));
    
    figure(38);clf; %set(gcf,'position',[791 37 883 944])
    subplot(211);
    if ~isempty(wd_mask),
        imagesc(theta/pi*180,Z(1:end/2,1),double(sqrt(real(n2in(1:end/2,isort))))); 
        hold on; axis ij
        plot(theta/pi*180,wd_mask(1,isort),'w','linewidth',2);
        ylim([ 0 max(wd_mask(1,:))*1.2]);
    end
    xlabel('Y (Deg)'); ylabel('Z (m)'); title(dista); colorbar
    set(gca,'tickdir','out')
    
    subplot(212)
    imagesc(theta/pi*180,Z(1:end/2,1),double(20*log10(abs(psi(1:end/2,isort))))); hold on; axis ij
    if ~isempty(wd_mask), 
        plot(theta/pi*180,wd_mask(1,isort),'w','linewidth',2);
        ylim([ 0 max(wd_mask(1,:))*1.2]); 
    end
    caxis([-60 0]+max(caxis)); colorbar
    xlabel('Y (Deg)'); ylabel('Z (m)'); title(dista); 
    set(gca,'tickdir','out')
    drawnow
    
end

return


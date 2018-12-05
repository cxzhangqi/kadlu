function [lon,lat] = sub_transfer_XY_to_LL(x,y,lon_ref,lat_ref,rot,x_ref,y_ref)

deg2rad = pi/180;

if nargin == 5,
    x_ref = 0;
    y_ref = 0;
elseif nargin == 4,
    x_ref = 0;
    y_ref = 0;
    rot = 0;  % counterclock angle
end

%----------------------------------
%  invert the planar approximation
%  x -- easten distance PLUS coord rotation
%  y -- north distance PLUS coord rotation
%----------------------------------

R = 6371009;        % Equatorial radius (6,378.1370 km)
                    % Polar radius (6,356.7523 km)
                    % The International Union of Geodesy and Geophysics 
                    %    (IUGG) defines the mean radius (denoted R1) to be                        
                    %    6,371.009 km
R2 = R*cos(lat_ref*deg2rad);

x = x-x_ref;
y = y-y_ref;

if rot ~= 0,
    SINROT = sin(rot*deg2rad);
    COSROT = cos(rot*deg2rad);
    XY = [COSROT -SINROT; SINROT COSROT] * [x(:).';   y(:).'];
    x(:) = XY(1,:);
    y(:) = XY(2,:);
end

lon = lon_ref + x/deg2rad/R2;   %X
lat = lat_ref + y/deg2rad/R;  %Y

return

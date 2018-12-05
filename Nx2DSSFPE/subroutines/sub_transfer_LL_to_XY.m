function [x,y] = sub_transfer_LL_to_XY(lon,lat,lon_ref,lat_ref,rot,x_ref,y_ref)

% [x,y] = sub_transfer_LL_to_XY(lon,lat,lon_ref,lat_ref,rot,x_ref,y_ref)

deg2rad = pi/180;

if nargin == 5,
    x_ref = 0;
    y_ref = 0;
elseif nargin == 4,
    x_ref = 0;
    y_ref = 0;
    rot = 0;  % counterclock angle
end

%------------------------
%  planar approximation
%  x -- easten distance PLUS coord rotation
%  y -- north distance PLUS coord rotation
%------------------------

R = 6371009;        % Equatorial radius (6,378.1370 km)
                    % Polar radius (6,356.7523 km)
                    % The International Union of Geodesy and Geophysics 
                    %    (IUGG) defines the mean radius (denoted R1) to be                        
                    %    6,371.009 km
R2 = R*cos(lat_ref*deg2rad);

x =(lon-lon_ref)*deg2rad*R2;   %X
y =(lat-lat_ref)*deg2rad*R;  %Y

if rot ~= 0,
    SINROT = sin(rot*deg2rad);
    COSROT = cos(rot*deg2rad);
    XY = [COSROT SINROT; -SINROT COSROT] * [x(:).';   y(:).'];
    x(:) = XY(1,:);
    y(:) = XY(2,:);
end

x = x+x_ref;
y = y+y_ref;

return


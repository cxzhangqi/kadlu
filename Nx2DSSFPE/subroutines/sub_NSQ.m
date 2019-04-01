function [f,dx,dy,dz] = sub_NSQ(x,y,z)

global ENV

%  Arbitrary 3-D sound speed field with trilinear interpolation
%
%  x, y, z, f, dx, dy, dz are all column vectors and in the same length!! 
%     and z positive pointing down
% 
%  ENV.NSQ.field = [x,y,z]
%  ENV.NSQ.dx, ENV.NSQ.x0, ENV.NSQ.x1
%  ENV.NSQ.dy, ENV.NSQ.y0, ENV.NSQ.y1
%  ENV.NSQ.dz, ENV.NSQ.z0, ENV.NSQ.z1
%
%  Y.-T. Lin 03/10/2014 @ WHOI




ENV.NSQ.field = ones(size(ENV.NSQ.field));




x = x(:); y = y(:); z = z(:);  % making sure every vector is in columns 
z = abs(z);  % z positive pointing down
tmp0 = zeros(size(x)); tmp1 = ones(size(x));  % zeros and ones

[nx,ny,nz] = size(ENV.NSQ.field); % check dimensions
if nx > 1,
    x(x<ENV.NSQ.x0)=ENV.NSQ.x0; x(x>=ENV.NSQ.x1)=ENV.NSQ.x1-1e-6;
    x = (x-ENV.NSQ.x0)/ENV.NSQ.dx;
    i = floor(x)+1; iplus1 = i+1; x = x - i +1;
else
    i = tmp1; iplus1 = tmp1; x = tmp0;
end
if ny > 1,
    y(y<ENV.NSQ.y0)=ENV.NSQ.y0; y(y>=ENV.NSQ.y1)=ENV.NSQ.y1-1e-6;
    y = (y-ENV.NSQ.y0)/ENV.NSQ.dy;
    j = floor(y)+1; jplus1 = j+1; y = y - j +1;
else
    j = tmp1; jplus1 = tmp1; y = tmp0;
end
if nz > 1,
    z(z<ENV.NSQ.z0)=ENV.NSQ.z0; z(z>=ENV.NSQ.z1)=ENV.NSQ.z1-1e-6;
    z = (z-ENV.NSQ.z0)/ENV.NSQ.dz;
    k = floor(z)+1; kplus1 = k+1; z = z - k +1;
else
    k = tmp1; kplus1 = tmp1; z = tmp0;
end

p000=squeeze(ENV.NSQ.field(sub2ind([nx ny nz],i,j,k)));
p001=squeeze(ENV.NSQ.field(sub2ind([nx ny nz],i,j,kplus1)));
p010=squeeze(ENV.NSQ.field(sub2ind([nx ny nz],i,jplus1,k)));
p011=squeeze(ENV.NSQ.field(sub2ind([nx ny nz],i,jplus1,kplus1)));
p100=squeeze(ENV.NSQ.field(sub2ind([nx ny nz],iplus1,j,k)));
p101=squeeze(ENV.NSQ.field(sub2ind([nx ny nz],iplus1,j,kplus1)));
p110=squeeze(ENV.NSQ.field(sub2ind([nx ny nz],iplus1,jplus1,k)));
p111=squeeze(ENV.NSQ.field(sub2ind([nx ny nz],iplus1,jplus1,kplus1)));

B1 = [1  0  0  0  0  0  0 0
     -1  0  0  0  1  0  0 0
     -1  0  1  0  0  0  0 0
     -1  1  0  0  0  0  0 0
      1  0 -1  0 -1  0  1 0
      1 -1 -1  1  0  0  0 0
      1 -1  0  0 -1  1  0 0
     -1  1  1 -1  1 -1 -1 1] ;

% C = B1*[p000 p001 p010 p011 p100 p101 p110 p111].'
% Q1 = [1 x y z xy yz zx xyz].'
% p = C'*Q1 = [p000 p001 p010 p011 p100 p101 p110 p111]*B1'*Q1

if nargout > 1, cal_grad = 1; else cal_grad = 0; end
 
% f = nan(size(x));
% if cal_grad,
%     dx = nan(size(x)); dy = nan(size(x)); dz = nan(size(x));
% end
% parfor ind = 1:length(x)
%     f(ind) = [p000(ind) p001(ind) p010(ind) p011(ind) p100(ind) p101(ind) p110(ind) p111(ind)]*B1.' ...
%             * [1 x(ind) y(ind) z(ind) x(ind)*y(ind) y(ind)*z(ind) z(ind)*x(ind) x(ind)*y(ind)*z(ind)].';
%     if cal_grad,
%         dx(ind) = [p000(ind) p001(ind) p010(ind) p011(ind) p100(ind) p101(ind) p110(ind) p111(ind)]*B1.' ...
%             * [0 1 0 0 y(ind) 0      z(ind) y(ind)*z(ind)].';
%         dy(ind) = [p000(ind) p001(ind) p010(ind) p011(ind) p100(ind) p101(ind) p110(ind) p111(ind)]*B1.' ...
%             * [0 0 1 0 x(ind) z(ind) 0      x(ind)*z(ind)].';
%         dz(ind) = [p000(ind) p001(ind) p010(ind) p011(ind) p100(ind) p101(ind) p110(ind) p111(ind)]*B1.' ...
%             * [0 0 0 1 0      y(ind) x(ind) x(ind)*y(ind)].';
%     end
% end
% if cal_grad,
%     if nx > 1, dx = dx./ENV.NSQ.dx; end
%     if ny > 1, dy = dy./ENV.NSQ.dy; end
%     if nz > 1, dz = dz./ENV.NSQ.dz; end
% end

% vectorized computation
B1 = [p000 p001 p010 p011 p100 p101 p110 p111]*B1.'; 
if cal_grad, 
    xy = x.*y; yz = y.*z; zx = z.*x;
     f = sum( B1.*[tmp1 x y z xy yz zx xy.*z] , 2); 
    dx = sum( B1.*[tmp0 tmp1 tmp0 tmp0 y tmp0 z yz] , 2); if nx > 1, dx = dx./ENV.NSQ.dx; end
    dy = sum( B1.*[tmp0 tmp0 tmp1 tmp0 x z tmp0 zx] , 2); if ny > 1, dy = dy./ENV.NSQ.dy; end
    dz = sum( B1.*[tmp0 tmp0 tmp0 tmp1 tmp0 y x xy] , 2); if nz > 1, dz = dz./ENV.NSQ.dz; end
else
    xy = x.*y;
    f = sum( B1.*[tmp1 x y z xy y.*z z.*x xy.*z] , 2); 
end

return
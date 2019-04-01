function [f,dy,dyy,n] = sub_SeafloorDepth(x,y)

global ENV

%  Arbitrary Sea floor depth with bilinear interpolation
%  
%  x, y, f, dy, dyy, n are all column vectors and in the same length!! 
% 
%  ENV.WD.field = [x,y]  % positive water depth, negative landscape
%  ENV.WD.dx, ENV.WD.x0, ENV.WD.x1
%  ENV.WD.dy, ENV.WD.y0, ENV.WD.y1
%
%  Y.-T. Lin 03/10/2014 @ WHOI

x = x(:); y = y(:); % making sure every vector is in columns 
tmp0 = zeros(size(x)); tmp1 = ones(size(x));  % zeros and ones



ENV.WD.field = ones(size(ENV.WD.field)) * 10000;



[nx,ny] = size(ENV.WD.field); % check dimensions 
if nx > 1,
    x(x<ENV.WD.x0)=ENV.WD.x0; x(x>=ENV.WD.x1)=ENV.WD.x1-1e-6;
    x = (x-ENV.WD.x0)/ENV.WD.dx;
    i = floor(x)+1; iplus1 = i+1; x = x + 1 - i;
else
    i = tmp1; iplus1 = tmp1; x = tmp0;
end
if ny > 1;
    y(y<ENV.WD.y0)=ENV.WD.y0; y(y>=ENV.WD.y1)=ENV.WD.y1-1e-6;
    y = (y-ENV.WD.y0)/ENV.WD.dy;
    j = floor(y)+1; jplus1 = j+1; y = y + 1 - j;
else
    j = tmp1; jplus1 = tmp1; y = tmp0;
end

p00 = squeeze(ENV.WD.field(sub2ind([nx ny],i,j)));
p01 = squeeze(ENV.WD.field(sub2ind([nx ny],i,jplus1)));
p10 = squeeze(ENV.WD.field(sub2ind([nx ny],iplus1,j)));
p11 = squeeze(ENV.WD.field(sub2ind([nx ny],iplus1,jplus1)));

B1 = [ 1  0  0  0
      -1  0  1  0
      -1  1  0  0
       1 -1 -1  1 ] ;

% bilinear interpolation
% ----------------------
% p(x,y) = p00 + (p10 - p00)[(x - x0)/(x1 - x0)]
%              + (p01 - p00)[(y - y0)/(y1 - y0)]
%              + (p11 - p01 - p10 + p00)[(x - x0)/(x1 - x0)][(y - y0)/(y1 - y0)], 
%
% which can be determined from the following procedure. 
% 
% C = B1*[p00 p01 p10 p11].';  % which means C = [p00 p10-p00 p01-p00 p11-p01-p10+p00] 
% let Q1 = [1 x y xy].'
% p = C.'*Q1 = [p00 p01 p10 p11]*B1.'*Q1

if nargout > 1, cal_grad = 1; else cal_grad = 0; end

% f = nan(size(x));
% if cal_grad,
%     dx = nan(size(x)); dy = nan(size(x));
% end
% parfor ind = 1:length(x)
%     f(ind) = [p00(ind) p01(ind) p10(ind) p11(ind)]...
%         *B1.'...
%         *[1 x(ind) y(ind) x(ind)*y(ind)].';
%     if cal_grad,
%         dx(ind) = [p00(ind) p01(ind) p10(ind) p11(ind)]...
%             *B1.'...
%             *[0 1 0 y(ind)].';
%         dy(ind) = [p00(ind) p01(ind) p10(ind) p11(ind)]...
%             *B1.'...
%             *[0 0 1 x(ind)].';
%     end
% end
% if cal_grad,
%     if nx > 1, dx = dx./ENV.WD.dx; end
%     if ny > 1, dy = dy./ENV.WD.dy; end
%     dyy = 0*dy; 
% end

% vectorized computation
B1 = [p00 p01 p10 p11]*B1.'; 
f = sum( B1.*[tmp1 x y x.*y] , 2); 
if cal_grad, 
    dx = sum( B1.*[tmp0 tmp1 tmp0 y] , 2); if nx > 1, dx = dx./ENV.WD.dx; end
    dy = sum( B1.*[tmp0 tmp0 tmp1 x] , 2); if ny > 1, dy = dy./ENV.WD.dy; end
    dyy = 0*dy; 
end

if nargout == 4,
    dz = -1*ones(size(x));
    norm = sqrt(dx.^2+dy.^2+dz.^2);
    n = [dx./norm; dy./norm; dz./norm]; % normalization
end


return


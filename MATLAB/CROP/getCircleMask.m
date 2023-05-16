function mask = getCircleMask(I,center,radius)



thetaResolution = 2; 
theta=(0:thetaResolution:360)'*pi/180;

x = bsxfun(@times,radius',cos(theta));
x = bsxfun(@plus,x,(center(:,1))');
x = cat(1,x,nan(1,length(radius)));
x = x(:);
x = x(~isnan(x));

y = bsxfun(@times,radius',sin(theta));
y = bsxfun(@plus,y,(center(:,2))');
y = cat(1,y,nan(1,length(radius)));
y = y(:);
y = y(~isnan(y));

mask = roipoly(I,x,y);
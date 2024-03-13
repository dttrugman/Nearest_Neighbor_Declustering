function D = latlonkm(p1,p2);

% Surface distance (in km) between points given as [lat,lon,depth]
% If depth is not specified, the surface distance is computed
%
% Usage
%     D = latlonkm([lat1 lon1],[lat2 lon2]);
%     D = latlonkm([lat1 lon1 depth1],[lat2 lon2 depth2]);
%

lat1=p1(:,1);
lon1=p1(:,2);
lat2=p2(:,1);
lon2=p2(:,2);

if size(p1,2)==3
    dep1=p1(:,3)*1000; % km -> m
end
if size(p2,2)==3
    dep2=p2(:,3)*1000;
end

R0 = 6.3673e6; % Earthquake radius in meters
D = R0*acos(sin(lat1*pi/180.0).*sin(lat2*pi/180.0)...
         + cos(lat1*pi/180.0).*cos(lat2*pi/180.0).*cos((lon1 - lon2)*pi/180.0));
if exist('dep1')==1
     D = sqrt(D.^2+(dep1-dep2).^2);
end
D=D/1000; % m -> km
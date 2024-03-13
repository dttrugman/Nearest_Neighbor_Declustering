function [P,n]=bp_2cat_add_1(time0,mag0,Lon0,Lat0,depth0,time,mag,Lon,Lat,depth,b,df,lag,iswaitbar);

%% Nearest-neighbor Baiesi-Paczuski proximity 
% for events in catalog (time,Lon,Lat,mag,depth) with respect
% to events from catalog (time0,Lon0,Lat0,mag0,depth0)
% Faster algorithm (see also bp_add_1.m) 

%% Parameters
%========================================
tmin=1; % seconds
dmin=1; % meters
%========================================

if exist('iswaitbar')~=1
    iswaitbar = 1;
end

[tmp,imin]=min(time);
P(imin,1)=0; % first element is the root
n(imin,1)=Inf;
T(imin,1)=NaN;
D(imin,1)=NaN;
M(imin,1)=NaN;

y2s = 365.25*24*60*60; % years to seconds
time=time*y2s;
time0=time0*y2s;
if exist('lag')~=1
    lag=Inf;
else
    lag=lag*y2s;
end

kl=cos(mean(Lat0)/180*pi);

if iswaitbar
    wbh=waitbar(0,'Please wait...');
    set(wbh,'Name','Baiesi-Paczuski metric on 2 catalogs, 2018')
end
di=0;
dis=0;
for i=1:length(time)
    %I=find(time0<time(i) & time0>=time(i)-lag & ...
    %    ~(Lat0==Lat(i) & Lon0==Lon(i)));
    
    % Select a small number of potential parents
    %------------------------------------------
    % Upper estimate for d to previous event
    ind = max(find(time0<time(i))); % immediate previous event in cat0
    K = max(1,ind-50):ind; % 50 earlier events in cat0
    if ~isempty(K)
        kl0 = cos(mean(Lat0(K))/180*pi); % latitudinal corection
        d = sqrt(((Lon(i)-Lon0(K))*kl0).^2+(Lat(i)-Lat0(K)).^2)*111; % distance to i in cat
        d = sqrt(d.^2 + (depth(i)-depth0(K)).^2);
        d = max(d*1000,dmin)/1000; % in km

        % Time difference
        t = max(time(i)-time0(K),tmin)/y2s; % in years
        
        
        J = find(Lon(i)~=Lon0(K) | Lat(i)~=Lat0(K) | depth(i)~=depth0(K));
        t = t(J);
        d = d(J);
        K = K(J);
        
        % log-proximity to previous event
        nc0 = min(log10(t)+df*log10(d)-b*mag0(K)); % years, km
    
        % Point "in the middle" of nc0 constant-proximity level
        q = 3.36; % empirical, can be estimated by a single run of bp_add_1.m 
        dt0 = 10.^((nc0-q)/2);
        d0 = 10.^((nc0+q)/2/df);
        
        % All events from cat0 within dt0 or d0
        % Coefficient 2 is > sqrt(3) = 1.73 needed in 3D space
        I = find((time(i)-time0(1:ind)   <= dt0*y2s | ...         % within dt0
            (abs(Lon(i)-Lon0(1:ind))     <= 2*d0/111/kl & ...  % within d0 in Lon
             abs(Lat(i)-Lat0(1:ind))     <= 2*d0/111 & ...   % within d0 in Lat
             abs(depth(i)-depth0(1:ind)) <= 2*d0)) & ...         % withid d0 in depth
            ~(Lon0(1:ind)==Lon(i) & Lat0(1:ind)==Lat(i) & depth0(1:ind)==depth(i))); % not the same as i from cat
    else
        I=[];
    end

    if ~isempty(I)
        if 1
            %d=latlonkm([Lat(i) Lon(i)],[Lat0(I) Lon0(I)]);
            d=latlonkm([Lat(i) Lon(i) depth(i)],[Lat0(I) Lon0(I) depth0(I)]);
            if di==0
                 disp('Geographic distance');di=1;
            end
            
        end  
 
        if 0
            %d=sqrt((Lat(i)-Lat0(I)).^2+(Lon(i)-Lon0(I)).^2);
            d=sqrt((Lat(i)-Lat0(I)).^2+(Lon(i)-Lon0(I)).^2+(depth(i)-depth0(I)).^2);
            if di==0
                disp('Euclidean distance');di=1;
            end
            
        end
        
        d=d*1000; % meters
        d=max(d,dmin)/1000; % km
        t=max(time(i)-time0(I),tmin)/365.25/24/60/60; % years
        nc = log10(t)+df*log10(d)-b*mag0(I);
        
        [~,imin]=min(nc);

        n(i,1)=nc(imin);
        P(i,1)=I(imin);
        D(i,1)=d(imin);
        T(i,1)=t(imin);
        M(i,1)=mag0(I(imin));
    else
        n(i,1)=Inf;
        P(i,1)=0;
        D(i,1)=0;
        T(i,1)=0;
        M(i,1)=0;
    end
    if iswaitbar
        if mod(i,10)==0
            waitbar(i/length(time));
        end
    end
end

if iswaitbar
    close(wbh);
end

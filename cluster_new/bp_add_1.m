function [P,n,D,T,M,n1,nc1]=bp_add_1(time,mag,Lon,Lat,depth,b,df,lag,iswaitbar)

%% Baiesi-Paczuski metric [PHYSICAL REVIEW E 69, 066106 (2004)]
% Additive verion: calculates log10(eta), which should run faster
%
% Inputs: 
%       time is the ocurrence time of events in years
%       mag  is the event magnitude
%       (Lon,Lat,depth) are event coordinates
%       b    is the GR slope
%       df   is the fractal dimension of epicenters
% Outputs (all are the same size as the EQ catalog):
%      P is the resulting tree (parent-pointer format)
%      n is the nearest-neighbor BP-distance
%      D is the nearest-neighbor spatial distance
%      T is the nearest-neighbor temporal distance
%      M is the magniude of the parent event
%
% January 22, 2021: distance lag

%% Parameters
%========================================
tmin=1; % seconds
dmin=1; % meters
%========================================

%% Sort by time
[time,Is] = sort(time);
mag = mag(Is);
Lon = Lon(Is);
Lat = Lat(Is);
if ~isempty(depth)
    depth = depth(Is);
end

%% Set up
if exist('iswaitbar')~=1
    iswaitbar = 1;
end

P(1,1)=0; % first element is the root
n(1,1)=Inf;
T(1,1)=NaN;
D(1,1)=NaN;
M(1,1)=NaN;

%% Convert time to seconds, set up lag
y2s = 365.25*24*60*60; % years to seconds
time=time*y2s; % yr -> sec
if exist('lag')~=1
    lag=Inf;
else
    lag=lag*y2s; % in sec
end

kl=cos(mean(Lat)/180*pi);

%% Main loop
if iswaitbar
    wbh=waitbar(0,'Please wait...');
    set(wbh,'Name','Baiesi-Paczuski proximity, 2004, 2021')
end
di=0;
dis=0;
for i = 2:length(time)
    if i<1000
         nc0 = Inf;
         I=find(time<time(i) & time>=time(i)-lag);
    else
        % Select a small number of potential parents
        %------------------------------------------
        % Upper estimate for d to previous event
        K = i+[-31:-1]; % 30 earlier events
        kl0 = cos(mean(Lat([i K]))/180*pi); 
        d = sqrt(((Lon(i)-Lon(K))*kl0).^2+(Lat(i)-Lat(K)).^2)*111; % km
        if ~isempty(depth)
            d = sqrt(d.^2 + (depth(i)-depth(K)).^2); % km
        end
        d = max(d*1000,dmin)/1000; % km

        % Time difference
        t = max(time(i)-time(K),tmin)/y2s; % years

        % log-proximity to previous event
        nc0 = min(log10(t)+df*log10(d)-b*mag(K)); % years, km

        q = 3.36;
        dt0 = 10.^((nc0-q)/2);
        d0 = 10.^((nc0+q)/2/df);

        I = find(time(i)-time(1:i-1)    <= dt0*y2s |...
             abs(Lon(i)-Lon(1:i-1))     <= 1.1*d0/111/kl & ...
             abs(Lat(i)-Lat(1:i-1))     <= 1.1*d0/111 & ...
             abs(depth(i)-depth(1:i-1)) <= 1.1*d0);

        n1(i) = length(I);
        nc1(i) = nc0;
    end

    if ~isempty(I)
        if 1
            % possibility of events at the same locations
            J = find(Lat(i)~=Lat(I) | Lon(i)~=Lon(I) | depth(i)~=depth(I)); % not same locations
            d = zeros(size(Lat(I)));
            if ~isempty(J)
                d(J)=latlonkm([Lat(i) Lon(i) depth(i)],[Lat(I(J)) Lon(I(J)) depth(I(J))]);
            end
            if di==0
                disp('Geographic distance');di=1;
            end
        end
        
        if 0
            %d=sqrt((Lat(i)-Lat(I)).^2+(Lon(i)-Lon(I)).^2);
            d=sqrt((Lat(i)-Lat(I)).^2+(Lon(i)-Lon(I)).^2+(depth(i)-depth(I)).^2);
            if di==0
                disp('Euclidean distance');di=1;
            end
        end   
        d=d*1000; %meters
        
        d=max(d,dmin)/1000; % km
        t=max(time(i)-time(I),tmin)/y2s; % years
        nc=log10(t)+df*log10(d)-b*mag(I);

        [n(i,1),imin]=min(nc);
        P(i,1)=I(imin);
        D(i,1)=d(imin);
        T(i,1)=t(imin);
        M(i,1)=mag(I(imin));
    else
        n(i,1)=Inf;
        P(i,1)=0;
        D(i,1)=NaN;
        T(i,1)=NaN;
        M(i,1)=NaN;
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

function [a,eta,Lbm,P,D,T,M,Lb,Q] = bp_thinning_fast(time,mag,Lon,Lat,depth,b,df,eta0,Nboot,lag)

%% Reshuffled threshold for earthquake distance 
% Version October 30, 2021 

Ikeep=ones(size(time));

%% Original catalog
%---------------------------------------
[P,eta,D,T,M]=bp_add_1(time,mag,Lon,Lat,depth,b,df,lag,1);
J = find(eta(:,1)>eta0);

%% Permuted catalog
%---------------------------------------
Q=zeros(size(time)); % bootstrap quantiles for eta 
Lb=zeros(size(time)); % linear bootstrap eta
Le=zeros(size(time)); % logarithmic bootstrap eta
M2b=zeros(size(time)); % 2nd moment, for variance calculation
Nboot_actual = zeros(size(time));

wbh=waitbar(0,'Reshuffled catalog analysis. Please wait...');
set(wbh,'Name','Earthquake declustering: nearest neighbor proximity analysis (c) 2018');
for iboot = 1:Nboot
    
    e_dur = 0;
    % Extend catalog before the first event (to treat correctly the initial part)
    % for duration e_dur
    if e_dur>0
        Nc = length(J); % current numer of treated events
        dur = range(time) + e_dur; % add 5 years
        N = round(Nc*dur/range(time)); % new total number of events
        Next = N-Nc; % number of new events
        range_ext = [min(time)-e_dur max(time)];  % new time interval
        Iext = randperm(Nc);
        Iext = Iext(1:Next); % random sample of length Next from J
        Lon_ext = Lon(J(Iext)); % new Longitudes
        Lat_ext = Lat(J(Iext)); % new Latitudes
        Iext = randperm(Nc);
        Iext = Iext(1:Next); % random sample of length Next from J
        mag_ext = mag(J(Iext)); % new magnitudes

        Lons = [Lon(J);Lon_ext];
        Lats = [Lat(J);Lat_ext];
        mags = [mag(J);mag_ext];
    else
        Lons = Lon(J);
        Lats = Lat(J);
        mags = mag(J);
        N = length(J);
        range_ext = [min(time) max(time)];
    end
    
    % Random/Reshuffled catalog parameters 
    if 0
        % Reshuffled times
        timeb = time(J(randperm(length(J)))); % reshuffled time
    elseif 1
        % Uniform time
        timeb = sort(rand(N,1)*range(range_ext)+min(range_ext)); % uniform time
    elseif 0
        % Poisson process with kernel intensity
        [ma,nn,ss,wa,v,maa,nw,tout] = mov_ave(time(J),mag(J),10,time);
        timeb = ppi_conditioned(nn',tout',N);
    end
    % random jitter is necessary to prevent 0-distance from events of original catalog  
    Ip = randperm(N);
    Lonb = Lons(Ip);
    Latb = Lats(Ip);
    depthb = depth(J(Ip));
    magb = mags(randperm(N));
    [Pb,etab]=bp_2cat_add_1(timeb,magb,Lonb,Latb,depthb,time,mag,Lon,Lat,depth,b,df,lag,0);
    
    I=find(etab < Inf);
    Lb(I) = Lb(I) + etab(I); % logarithmic version
    M2b(I) = M2b(I) + etab(I).^2; % 2nd moment, for variance calculation 
    Nboot_actual(I) = Nboot_actual(I)+1;
    
    Q = Q + 1*(etab<eta);

    waitbar(iboot/Nboot);
end
close(wbh);

Lbm = Lb./Nboot_actual;
Sdb=ones(size(M2b));
Q = Q/Nboot;
a = (eta-Lbm)./Sdb;

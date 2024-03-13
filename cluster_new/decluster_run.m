%% Declustering test script
% Run it to test that all dependencies are present

%% Parameters
df = 2.6;    % fractal dimension of hypocenters 
eta0 = 1;   % initial threshold to remove most clustered events
Nboot = 100; % number of catalog reshufflings
alpha0 = 0;  % cluster threshold
%=============

clear
load input_catalog
[ad0,etad0,Lbmd0,Pd0,Dd0,Td0,Md0,Lbd0,Qd0] = bp_thinning_fast(time,mag,Lon,Lat,depth,0,df,eta0,Nboot,Inf);
cluster_analysis
save declustered_catalog

%% Plot results
figure(1)
subplot(2,1,1)
plot(time,Lat,'.','markersize',1,'markerfacecolor','b');
%scatter(time,Lat,10.^(0.75*(mag-min(mag)+.5)),'o','filled','markeredgecolor','k')
xlabel('Time, years')
ylabel('Latitude')
title('All events')
allx([min(time) max(time)]);
setfig

subplot(2,1,2)
plot(time(Imain),Lat(Imain),'.','markersize',1,'markerfacecolor','b');
%scatter(time(Imain),Lat(Imain),10.^(0.75*(mag(Imain)-min(mag)+.5)),'o','filled','markeredgecolor','k')
xlabel('Time, years')
ylabel('Latitude')
title('Background events')
allx([min(time) max(time)]);
setfig
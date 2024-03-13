%% Declustering test script
% Run it to test that all dependencies are present

clear
load ANSS_test
df = 2.6;    % fractal dimension of hypocenters 
eta0 = 1;    % initial threshold to remove most clustered events
Nboot = 10;  % number of catalog reshufflings
alpha0 = 0.5;  % cluster threshold

[ad0,etad0,Lbmd0,Pd0,Dd0,Td0,Md0,Lbd0,Qd0] = bp_thinning_fast(time,mag,Lon,Lat,depth,0,df,eta0,10,Inf);
cluster_analysis

%% Plot results
figure(1)
subplot(2,1,1)
plot(time,Lat,'ko','markersize',4,'markerfacecolor','b');
%scatter(time,Lat,10.^(0.75*(mag-min(mag)+.5)),'o','filled','markeredgecolor','k')
xlabel('Time, years')
ylabel('Latitude')
title('All events')
allx([min(time) max(time)]);
setfig

subplot(2,1,2)
plot(time(Imain),Lat(Imain),'ko','markersize',4,'markerfacecolor','b');
%scatter(time(Imain),Lat(Imain),10.^(0.75*(mag(Imain)-min(mag)+.5)),'o','filled','markeredgecolor','k')
xlabel('Time, years')
ylabel('Latitude')
title('Background events')
allx([min(time) max(time)]);
setfig
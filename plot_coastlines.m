
load('transects_with_shorelines_final.mat')
% TODO - automate lookup into list of littoral cells.
% For now: 
% index into all Core Banks transects
idx = 23820:25052;

% index into N. Core Banks only
idx = 24476:25052;
lc_name = 'Portsmouth Island'

% % Rodanthe (probably need to narrow this down)
% idx= 25896:27146;
% lc_name = 'Rodanthe'
%% Colors from ColorBrewer
load_colors

%% Calculate some shoreline location statistics
% declare an array for statistics
% stats = [mean, std, 25th, 50th, 75th, slope, std.err.slope, r2]
stats = nan*ones(length(idx),8);
utmxstats = nan*ones(length(idx),4);
utmystats = nan*ones(length(idx),4);
n = nan*ones(length(idx));
sdate = nan*ones(length(idx));
edate = nan*ones(length(idx));

% loop through indices for this location and calc stats
% TODO - adjust to look at specific time intervals
for i=1:length(idx)
    jdx = find([transects(idx(i)).SAT(:)]==1); % satellite shorelines

    Y = [transects(idx(i)).Y(jdx)];
    t = [transects(idx(i)).t(jdx)];
    n(i) = length(jdx);
    sdate(i) = min(t);
    edate(i) = max(t);
    % Calculate shoreline location stats in transect coordinates
    stats(i,1) = mean(Y);
    stats(i,2) = std(Y);
    tmp = prctile(Y,[25,50,75]); % needs statistics toolbox
    stats(i,3)=tmp(1);
    stats(i,4)=tmp(2);
    stats(i,5)=tmp(3);
    % linear fit to shoreline change
    [a,b,r2,sa,sb,hdot]=lsfit(t,Y,0); % replace with "official" USGS fitting method
    stats(i,6) = b*365.25; % convert slope to m/year
    stats(i,7) = sb; % std. error of slope
    stats(i,8) = r2;
    
    % Convert shoreline location stats to UTM
    % mean shoreline
    [utmxstats(i,1), utmystats(i,1)] = ...
        transect2utm( transects(idx(i)).x_on, transects(idx(i)).y_on, ...
        stats(i,1), transects(idx(i)).angle );
   % 25th pct shoreline
   [utmxstats(i,2), utmystats(i,2)] = ...
        transect2utm( transects(idx(i)).x_on, transects(idx(i)).y_on, ...
        stats(i,3), transects(idx(i)).angle );
   % 50th pct shoreline
   [utmxstats(i,3), utmystats(i,3)] = ...
        transect2utm( transects(idx(i)).x_on, transects(idx(i)).y_on, ...
        stats(i,4), transects(idx(i)).angle );
   % 75th pct shoreline
   [utmxstats(i,4), utmystats(i,4)] = ...
        transect2utm( transects(idx(i)).x_on, transects(idx(i)).y_on, ...
        stats(i,5), transects(idx(i)).angle );
end
%% Calc stats on stats
% shoreline change rates
fprintf(1,'Percentiles (5,25,50,75,95) of shoreline change rates for %s\n',lc_name)
LTEpct=prctile([transects(idx).LTER],[5,25,50,75,95])
fprintf(1,'N = %d\n',length(idx))
% this will be slightly different bc it does not include the non-satellite
% shorelines, and the fitting procedure might not be right
LTE2pct=prctile(stats(:,6),[5,25,50,75,95])
%% Plot median and quartile shoreline locations in UTM
% (you will have to zoom in to see them)
figure(1); clf
plot(utmxstats(:,3),utmystats(:,3),'.','color',purples(4,:))
hold on
plot(utmxstats(:,2),utmystats(:,2),'.','color',purples(2,:))
plot(utmxstats(:,4),utmystats(:,4),'.','color',purples(2,:))
xlabel('UTM x (m)')
ylabel('UTM y (m)')
title(lc_name)
shg
%%
figure(2); clf
plot(idx, zeros(size(idx)),'--k')
hold on
scatter(idx,stats(:,6),14,stats(:,8),'filled')
xlabel('Transect Number')
ylabel('<- Erosion (m/y) Depostion ->')
h = colorbar;
set(get(h,'label'),'string','Regression coefficient {\itr}^2');
title(lc_name)



%% load the data
load('c:\crs\proj\2020_CALO\FloSup_coastlines\transects_with_shorelines_final.mat')
%% Pick a section
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

%% Make an cumulative alongshore distance array for plotting
utmxy_Y0 = nan*ones(length(idx),2);
% find location of Y0
for i = 1:length(idx)
    [utmxy_Y0(i,1), utmxy_Y0(i,2)] = ...
        transect2utm( transects(idx(i)).x_on, transects(idx(i)).y_on, ...
        transects(idx(i)).Y0, transects(idx(i)).angle );
end
dxy = diff(utmxy_Y0);
% alongshore distance in km
dalong = [0; cumsum(sqrt( dxy(:,1).^2 + dxy(:,2).^2 ))]/1000.;
%% Save the UTM coords of Y0
% These are in UTM zone 17N
fnutmY0 = sprintf('%s_Y0.csv',lc_name)
fid = fopen(fnutmY0,'w');
for i=1:length(idx)
    fprintf(fid,'%12.2f,%12.2f\n',utmxy_Y0(i,:));
end
fclose(fid);
%% write out Y0 start locations for bounding boxes
bbx_list = unique([transects(idx).coastsat_bbox_number])';
% and make an array of bbox start points in dalong coords
dalong_bbox = nan*ones(length(bbx_list),2);
for i =1:length(bbx_list)
    ibx = find([transects(idx).coastsat_bbox_number]==bbx_list(i),1,'first');
    dalong_bbox(i,:) = [bbx_list(i), dalong(ibx)];
    fbbx = sprintf('bbx%05d.csv',bbx_list(i))
    fid = fopen(fbbx,'w');
    fprintf(fid,'%12.2f,%12.2f\n',utmxy_Y0(ibx,:))
    fclose(fid);
end
%% Calculate some shoreline location statistics for entire time
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
    % t statistic for alpha=0.5, two-tailed, DOF = n(i)-2
    tstat = tinv(1-.05/2,n(i)-2.)
    stats(i,7) = tstat*sb*365.25; % 90% conf. limits on slope
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


%% Calc shoreline change for specific time intervals
% Define time intervals
% Looking at 1984 to 2021 inclusive
dt = 13 % length of time to fit shoreline change rates
nint = ceil((2021-1984)/dt)
ystart = (1984:dt:1984+nint*dt)';
dnstart = datenum(ystart,ones(size(ystart)),ones(size(ystart)));
datestr(dnstart)

% declare arrays for statistics
% stats = [mean, std, 25th, 50th, 75th, slope, std.err.slope, r2]
stats_t = nan*ones(nint,length(idx),8);
utmxstats_t = nan*ones(nint,length(idx),4);
utmystats_t = nan*ones(nint,length(idx),4);
n_t = nan*ones(nint,length(idx),1);
sdate_t = nan*ones(nint,length(idx),1);
edate_t = nan*ones(nint,length(idx),1);

LTE2pct_t = nan*ones(nint,5);

% loop through time intervals
for k = 1:length(dnstart)-1
    fprintf(1,'%s to %s\n',datestr(dnstart(k)),datestr(dnstart(k+1)))

    % loop through transects and calc stats
    for i=1:length(idx)
        %jdx = find([transects(idx(i)).SAT(:)]==1); % satellite shorelines
        % TODO - This does not exclude non-satellite data
        jdx = [find([transects(idx(i)).t(:)]>=dnstart(k),  1,'first'):...
               find([transects(idx(i)).t(:)]< dnstart(k+1),1,'last')]; %time interval
        N = length(jdx);
        %fprintf(1,'%s to %s: N=%d\n',datestr(dnstart(k)),datestr(dnstart(k+1)),n)
        if(N>=2)
            Y = [transects(idx(i)).Y(jdx)];
            t = [transects(idx(i)).t(jdx)];
            n_t(k,i) = N;
            sdate_t(k,i) = min(t);
            edate_t(k,i) = max(t);
            % Calculate shoreline location stats in transect coordinates
            stats_t(k,i,1) = mean(Y);
            stats_t(k,i,2) = std(Y);
            tmp = prctile(Y,[25,50,75]); % needs statistics toolbox
            stats_t(k,i,3)=tmp(1);
            stats_t(k,i,4)=tmp(2);
            stats_t(k,i,5)=tmp(3);
            % linear fit to shoreline change
            [a,b,r2,sa,sb,hdot]=lsfit(t,Y,0); % replace with "official" USGS fitting method
            % t statistic for alpha=0.5, two-tailed, DOF = n(i)-2
            % (needs statistics toolbox)
            tstat = tinv(1-.05/2,n_t(k,i)-2.);
            stats_t(k,i,6) = b*365.25; % convert slope to m/year
            stats_t(k,i,7) = tstat*sb*365.25; % 90% confidence interval on slope
            stats_t(k,i,8) = r2;

            % Convert shoreline location stats to UTM
            % mean shoreline
            [utmxstats_t(k,i,1), utmystats_t(k,i,1)] = ...
                transect2utm( transects(idx(i)).x_on, transects(idx(i)).y_on, ...
                stats_t(k,i,1), transects(idx(i)).angle );
            % 25th pct shoreline
            [utmxstats_t(k,i,2), utmystats_t(k,i,2)] = ...
                transect2utm( transects(idx(i)).x_on, transects(idx(i)).y_on, ...
                stats_t(k,i,3), transects(idx(i)).angle );
            % 50th pct shoreline
            [utmxstats_t(k,i,3), utmystats_t(k,i,3)] = ...
                transect2utm( transects(idx(i)).x_on, transects(idx(i)).y_on, ...
                stats_t(k,i,4), transects(idx(i)).angle );
            % 75th pct shoreline
            [utmxstats_t(k,i,4), utmystats_t(k,i,4)] = ...
                transect2utm( transects(idx(i)).x_on, transects(idx(i)).y_on, ...
                stats_t(k,i,5), transects(idx(i)).angle );
        else
            n_t(k,i)=N;
            stats_t(k,i,:)=nan;
            utmxstats_t(k,i,:)=nan;
            utmystats_t(k,i,:)=nan;
        end
    end
    % Calc stats on stats
    % shoreline change rates
    fprintf(1,'N = %d\n',n_t(i))
    % this will be slightly different bc it does not include the non-satellite
    % shorelines, and the fitting procedure might not be right
    LTE2pct_t(k,:)=prctile(stats_t(k,:,6),[5,25,50,75,95])
end
%% Plot median and quartile shoreline locations in UTM
% (you will have to zoom in to see them)
figure(1); clf
plot(utmxstats(:,3),utmystats(:,3),'.','color',purples(4,:))
hold on
plot(utmxy_Y0(:,1),utmxy_Y0(:,2),'.','color',blues(3,:))
plot(utmxstats(:,2),utmystats(:,2),'.','color',purples(2,:))
plot(utmxstats(:,4),utmystats(:,4),'.','color',purples(2,:))
xlabel('UTM x (m)')
ylabel('UTM y (m)')
title(lc_name)
shg
%% Plot LTER (from satellites only)
figure(2); clf
plot(dalong, zeros(size(idx)),'--k')
hold on
plot(dalong_bbox(:,2),zeros(size(dalong_bbox(:,2))),'xr')

plot(dalong,stats(:,6)+stats(:,7),'-','color',purples(3,:))
plot(dalong,stats(:,6)-stats(:,7),'-','color',purples(3,:))

scatter(dalong,stats(:,6),14,stats(:,8),'filled')
ylim([-6,6])
xlabel('Alongshore distance (km)')
ylabel('<- Erosion (m/y) Depostion ->')
h = colorbar;
set(get(h,'label'),'string','Regression coefficient {\itr}^2');
title(lc_name)
%% Shoreline change rate for various time intervals
ts = datestr(dnstart)
figure(3); clf
plot(dalong, zeros(size(idx)),'--k')
hold on
plot(dalong_bbox(:,2),zeros(size(dalong_bbox(:,2))),'xr')
for k=1:length(dnstart)-1
    hk(k) = plot(dalong,stats_t(k,:,6),'-','linewidth',2,'color',purples(k+1,:));
end
ylim([-10,10])
xlabel('Alongshore distance (km)')
ylabel('<- Erosion (m/y) Depostion ->')
tss = sprintf('%s, dt=%d y',lc_name,dt)
title(tss)
legend(hk,ts(1:nint,:),'location','northwest')
shg

%% Plot time series of 13-year fits
% series of shorelines and fit for a single transect
i = 280 % pick a transect number in idx
figure(4); clf
% loop through time intervals
for k = 1:length(dnstart)-1
    jdx = [find([transects(idx(i)).t(:)]>=dnstart(k),  1,'first'):...
        find([transects(idx(i)).t(:)]< dnstart(k+1),1,'last')]; %time interval
    N = length(jdx);

    Y = [transects(idx(i)).Y(jdx)];
    t = [transects(idx(i)).t(jdx)];
    % linear fit to shoreline change
    [a,b,r2,sa,sb,hdot]=lsfit(t,Y,1);
    hold on
end
jdx = find([transects(idx(i)).SAT(:)]==1); % satellite shorelines

Y = [transects(idx(i)).Y(jdx)];
t = [transects(idx(i)).t(jdx)];
% linear fit to shoreline change
[a,b,r2,sa,sb,hdot]=lsfit(t,Y,1);
datetick
ylabel('<- Erosion (m) Depostion ->')

%% Make a running nyr fit and plot
% there are some hard-wired parameters in this plot
i = 280 % pick a transect number in idx
ylen = [3, 7, 13];
ts=[' 3 years'; ' 7 years';'13 years']
cols = [grays(2,:);blues(3,:);purples(5,:)]
figure(5); clf
tz = [datenum(1985,1,1):datenum(2020,1,1)];
zz = zeros(size(tz))
plot(tz,zz,'--k')
hold on
for ik = 1:length(ylen)
    yfirst = 1984:2021-ylen(ik)
    nyr = nan*ones(length(yfirst),1);
    ymid = nyr;
    scr = nyr;
    scr90 = nyr;
    r2a = nyr;
    Na = nyr;
    for k = 1:length(yfirst)
        ystrt = datenum(yfirst(k),1,1);
        yend =  datenum(yfirst(k)+ylen(ik),1,1);
        ymid(k) = datenum(yfirst(k)+0.5*ylen(ik),1,1);

        jdx = [find([transects(idx(i)).t(:)]>=ystrt,  1,'first'):...
            find([transects(idx(i)).t(:)]< yend,1,'last')]; %time interval
        N = length(jdx);
        Y = [transects(idx(i)).Y(jdx)];
        t = [transects(idx(i)).t(jdx)];
        % linear fit to shoreline change
        [a,b,r2,sa,sb,hdot]=lsfit(t,Y,0);
        tstat = tinv(1-.05/2,N-2.);
        scr(k) = b*365.25; % convert slope to m/year
        scr90(k) = tstat*sb*365.25; % 90% confidence interval on slope
        r2a(k) = r2;
    end
    hk(ik)=plot(ymid, scr,'-','linewidth',3,'color',cols(ik,:))
    hold on
    plot(ymid, scr+scr90,'-','linewidth',1,'color',cols(ik,:));
    plot(ymid, scr-scr90,'-','linewidth',1,'color',cols(ik,:));
end
ylabel('<- Erosion (m/y) Depostion ->')
datetick
tss = sprintf('%s at alongshore = %4.1f km',lc_name,dalong(i))
title(tss)
legend(hk,ts(1:nint,:),'location','northwest')
shg


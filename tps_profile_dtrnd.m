function tps_profile_dtrnd

%%% 11/2020  Tohid Nozad Khalil
%%% plot LOS profile using stamp results
%%% this feature is added to ps_plot function of StaMPS
%%% this funtion is called by "profile" PushButton in left down corner of
%%% figure. Click in point you wish to start profile on double-click at
%%% the end point. aplly width of profile (in meter) in box below the 
%%% PushButton (default width is 100 m)

prompt = 'add width for your frofile in meter';
prof_width = input(prompt)*0.0005;			


hold on
load prof_dtrnd

%%% plot profile and create a rectangle that profile is located in the
%%% midlle of it
[p_lon,p_lat]=getpts;
plot(p_lon,p_lat,'r');
strt_lon=p_lon(1); strt_lat=p_lat(1);
end_lon=p_lon(2); end_lat=p_lat(2);

prof_width_ll = local2llh([0;prof_width],[strt_lon end_lon]);
prof_width_geo = sqrt((strt_lon-prof_width_ll(1))^2+(end_lon-prof_width_ll(2))^2);
theta = atand((end_lat-strt_lat)/(end_lon-strt_lon));
alpha = 90 - theta;
inc_lat = prof_width_geo*sind(alpha);
inc_lon = prof_width_geo*cosd(alpha);

polygon = [ strt_lon+inc_lon	strt_lat-inc_lat
            strt_lon-inc_lon	strt_lat+inc_lat
            end_lon-inc_lon     end_lat+inc_lat
            end_lon+inc_lon     end_lat-inc_lat ];
in = inpolygon(lonlat(:,1),lonlat(:,2),polygon(:,1),polygon(:,2));
index = find(in==1);

%%% extract lonlat and mean_v of profile
prof_ph=ph_disp(index);
lonlat_prof=lonlat(index,:);

%%% calculate projected distance of each point on x-axes 
prof_x=zeros(length(index),1);
for i = 1:length(index)
  dist=llh2local([lonlat_prof(i,1);lonlat_prof(i,2);0],[ strt_lon; strt_lat;0]);
  prof_x(i) = sqrt((dist(1)).^2 + dist(2).^2);
end

%%% plot data
figure
plot(prof_x,prof_ph,'o','LineWidth',0.5,'MarkerEdgeColor',[1 0 0], ...
    'MarkerFaceColor',[1 0 0], 'MarkerSize',3); 
pbaspect([10 4 1])
ylabel('LOS Velocity (mm/yr)')
xlabel('Distance (Km)')





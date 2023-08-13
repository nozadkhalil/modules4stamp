function tps_profile

%%% 11/2020  Tohid Nozad Khalil
%%% plot LOS profile using stamp results
%%% this feature is added to ps_plot function of StaMPS
%%% this funtion is called by "profile" PushButton in left down corner of
%%% figure. Click in point you wish to start profile on double-click at
%%% the end point. aplly width of profile (in meter) in box below the 
%%% PushButton (default width is 100 m)

prompt = 'add width for your profile in m';
prof_width = km2deg(input(prompt)/1000,'earth');


hold on
disp('please select which data you want to see profile;')
disp('meanV : mean velocity with no changes')
disp('dtrnd : mean velocity with removed trend')
disp('VS1   : Compare mean velocity with & without trend')
disp('VS2   : Compare mean velocity with & without trend/earthquake coseismic')
disp('grid  : mean velocity with new griding size')

prompt = 'Please mention a name (between quote marks);  ';
value_type = input(prompt);
vt = 0;
while vt < 1
  if strcmpi(value_type,'meanV')
      load prof
      vt = vt+1;
    elseif  strcmpi(value_type,'dtrnd')
      load ps_plot_dtrnd
      load ps2
      vt = vt+1;
    elseif  strcmpi(value_type,'VS1')
      mV  = load('prof');
      mVD = load('ps_plot_dtrnd');
      lonlat = mV.lonlat;
      vt = vt+1;
    elseif  strcmpi(value_type,'VS2')
      load psver
      psname=['./ps',num2str(psver)];
      phuwname=['./phuw',num2str(psver)];
      pheqname=['./rmEq'];
      sclaname=['./scla',num2str(psver)];

      ps     = load(psname);
      drop   = getparm('drop_ifg_index');
      lambda = getparm('lambda');
      unwrap_ifg_index=setdiff([1:ps.n_ifg],drop);
      day=ps.day(unwrap_ifg_index,:);
      lonlat = ps.lonlat;
      G=[ones(size(day)),day-ps.master_day]; 
      ref_ps=ps_setref;
      scla  = load(sclaname);
      if isfield(scla,'C_ps_uw')
         ph_scla = scla.ph_scla - repmat(scla.C_ps_uw,1,size(scla.ph_scla,2)); 
         else
         ph_scla = scla.ph_scla;
      end
      clear scla xy ij 
      mV  = tps_phDisp(phuwname,ph_scla,G,unwrap_ifg_index,ps.n_ps,ref_ps,lambda);
      mEq = tps_phDisp(pheqname,ph_scla,G,unwrap_ifg_index,ps.n_ps,ref_ps,lambda);
      mVD = load('ps_plot_dtrnd_rmEQ');
      clear ps ph_scla G unwrap_ifg_index ref_ps drop lambda *name psver
      vt = vt+1;
    elseif  strcmpi(value_type,'grid')
      load meanV_grd
      vt = vt+1;
    else
      vt = 0.1;
      disp('not recognize value type')
      prompt = 'Please type dataset name correctly (between quote marks);  ';
      value_type = input(prompt);
  end
end

fprintf('please select begin/end points for profile.\n') 
%%% plot profile and create a rectangle that profile is 
%%% located in the midlle of it
[p_lon,p_lat]=getpts;
fprintf('Start profile at => lon: %2.3f lat: %2.3f\n',p_lon(1),p_lat(1)) 
fprintf('End   profile at => lon: %2.3f lat: %2.3f\n',p_lon(2),p_lat(2)) 

prf = [p_lon(1) p_lat(1)
       p_lon(2) p_lat(2)];
plot(prf(:,1),prf(:,2),'r', 'LineWidth',1.5);
plot(prf(:,1),prf(:,2),'*r', 'MarkerSize',9);
[len,az]=distance('rh',prf(1,1),prf(1,2),prf(2,1),prf(2,2));
alpha = 90 + az;
rt = [sind(alpha)*prof_width/2 
      cosd(alpha)*prof_width/2];
prf_plgn = [prf(1,1)+rt(1)  prf(1,2)+rt(2)
            prf(1,1)-rt(1)  prf(1,2)-rt(2)
            prf(2,1)-rt(1)  prf(2,2)-rt(2)
            prf(2,1)+rt(1)  prf(2,2)+rt(2)
            prf(1,1)+rt(1)  prf(1,2)+rt(2)];
in = inpolygon(lonlat(:,1),lonlat(:,2),prf_plgn(:,1),prf_plgn(:,2));
lonlat_prf=lonlat(in,:);
%%% calculate projected distance of each point on x-axes 
prof_x=zeros(length(lonlat_prf),1);
for i = 1:length(lonlat_prf)
    dist=llh2local([lonlat_prf(i,1);lonlat_prf(i,2);0],[ prf(1,1); prf(1,2);0]);
    prof_x(i) = sqrt((dist(1)).^2 + dist(2).^2);
end

if ~strcmpi(value_type(1:2),'VS')
   %%% extract lonlat and mean_v of profile
   prof_ph=ph_disp(in);
   fID = fopen('prof.dat','w');
   for i=1:length(prof_x);
      fprintf(fID,'%f \t\t %f \t\t %f \t\t %f\n',...
      prof_x(i),lonlat_prf(i,1),lonlat_prf(i,2),prof_ph(i));
   end
   fclose(fID);
   %%% plot data
   figure
   plot(prof_x,prof_ph,...
   'o','MarkerEdgeColor','red','MarkerFaceColor','red')
   pbaspect([10 4 1])
   ylabel('LOS Velocity (mm/yr)')
   xlabel('Distance (Km)')
elseif strcmpi(value_type(1:2),'VS')
   if strcmpi(value_type,'VS2')
      meanV = mV(in);
      mVEq  = mEq(in);
      Dtrnd = mVD.ph_disp(in); 
      fID = fopen('prof_VS2.dat','w');
      for i=1:length(prof_x);
           fprintf(fID,'%f \t\t %f \t\t %f \t\t %f \t\t %f \t\t %f\n',...
           prof_x(i),lonlat_prf(i,1),lonlat_prf(i,2),meanV(i),Dtrnd(i),mVEq(i));
      end
      fclose(fID);
      figure
      hold on
      plot(prof_x,meanV,...
      'o','MarkerEdgeColor','red','MarkerSize',5,'MarkerFaceColor','red')
      plot(prof_x,Dtrnd,...
      'o','MarkerEdgeColor','blue','MarkerSize',3,'MarkerFaceColor','blue')
      plot(prof_x,mVEq,...
      'o','MarkerEdgeColor','black','MarkerSize',3,'MarkerFaceColor','black')
      pbaspect([10 4 1])
      ylabel('LOS Velocity (mm/yr)')
      xlabel('Distance (Km)')
      legend({'meanV','Dtrndt','EQremoved'},'Location','northwest')
      grid on
   elseif strcmpi(value_type,'VS1')
      meanV = mV.ph_disp(in);
      Dtrnd = mVD.ph_disp(in);  
      name=['prof_VS1.dat'];
      fID = fopen(name,'w');
      for i=1:length(prof_x);
         fprintf(fID,'%f \t\t %f \t\t %f \t\t %f \t\t %f\n',...
         prof_x(i),lonlat_prf(i,1),lonlat_prf(i,2),meanV(i),Dtrnd(i));
      end
      fclose(fID);
      figure
      hold on
      plot(prof_x,meanV,...
      'o','MarkerEdgeColor','red','MarkerSize',5,'MarkerFaceColor','red')
      plot(prof_x,Dtrnd,...
      'o','MarkerEdgeColor','blue','MarkerSize',5,'MarkerFaceColor','blue')
      pbaspect([10 4 1])
      ylabel('LOS Velocity (mm/yr)')
      xlabel('Distance (Km)')
      legend({'meanV','Dtrndt'},'Location','northwest')
      grid on
   end
end


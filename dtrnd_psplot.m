function [h_fig,lims,ifg_data_RMSE,h_axes_all]=dtrnd_psplot(value_type,varargin)

%   ======================================================================
%   11/2020 Tohid Nozadkhail: Add profile ploting feature
%	use the profile PushButton to call ps_rofile function
%	click in point you wish to start profile on double-click at the end point
%	aplly width of profile (in meter) in box below the PushButton (default width is 100 m)
%   12/2020 Tohid Nozadkhail: Add plot ATM corrected mean velocity not using phase value
%   use tps_plot('meanV-da','a_powerlaw')
%   01/2021 Tohid Nozadkhail: Add offset plot
%   02/2023 Tohid Nozadkhail: Add reference on GPS
%   ======================================================================


stdargin = nargin ; 
parseplotprm  

reference_flag=0;

if stdargin<1
    help ps_plot
    error('not enough input args')
end



load psver
psname=['./ps',num2str(psver)];
ps=load(psname);
day=ps.day;
master_day=ps.master_day;
xy=ps.xy;
lonlat=ps.lonlat;
n_ps=ps.n_ps;
n_ifg=ps.n_ifg;
master_ix=sum(day<master_day)+1;
ref_ps=0;    
clear psname


hold on
glii = strcmpi(value_type, {'gps', 'losgps', 'inFrmGps', 'inShpGps'});
if  sum(glii) > 0 || strcmpi(value_type,'IntSeis') 
    mn = min(lonlat);
    mx = max(lonlat);
    crnr = [mn;mx(1) mn(2) ;mx ;mn(1) mx(2)];
    hdng=getparm('heading');
    frm = rotateFrame(crnr, hdng); 
    frm = [ frm; frm(1,:)]; %%% frm is the coordinates of frame' corners
end

if  sum(glii) > 0
    Gdat = load('gps.dat');
    plot(frm(:,1), frm(:,2),'r')
    inFrmGPS = inpolygon(Gdat(:,1), Gdat(:,2), frm(:,1), frm(:,2));
    quiver(Gdat(inFrmGPS,1), Gdat(inFrmGPS,2), Gdat(inFrmGPS,3), Gdat(inFrmGPS,4))
    plot(Gdat(inFrmGPS,1), Gdat(inFrmGPS,2), '.r','markersize',10)
    pause(3)
end

if  strcmpi(value_type,'inShpGps') 
    disp('please select the area you want to model GPS on it .......')
    [plgnLon, plgnLat] = getpts;
    plgn = [plgnLon  plgnLat; plgnLon(1) plgnLat(1)];
    inShpGPS = inpolygon(Gdat(:,1), Gdat(:,2), plgn(:,1), plgn(:,2));
    plot(Gdat(inShpGPS,1), Gdat(inShpGPS,2), '.k','markersize',10)
    plot(plgn(:,1), plgn(:,2), 'r', 'linewidth', 1.5);
end

mrp = strcmpi(value_type, {'meanv', 'rmeq', 'phase'});
if sum(mrp) > 0 
    disp('please select reference area .......')
    [plgnLon, plgnLat] = getpts;
    plgn = [plgnLon  plgnLat; plgnLon(1) plgnLat(1)];
    plot(plgn(:,1), plgn(:,2), 'r', 'linewidth', 1.5)
end



plot_flag=1;
lims=[];
ref_ifg=0;
ifg_list=[];
n_x=0;
cbar_flag=0;
textsize=0;
textcolor=[0 0 0.004];
lon_rg=[];
lat_rg=[];
units=[];
%%% reference radius in case of external data to be plotted
ref_radius_data=1000;

n_y=0;

phuwname=['./phuw',num2str(psver)];
EQname=['./rmEq'];
phuwsbname=['./phuw_sb',num2str(psver)];
phuwsbresname=['./phuw_sb_res',num2str(psver)];
scnname=['./scn',num2str(psver)];
ifgstdname=['./ifgstd',num2str(psver)];
apsbandsname = ['./tca_bands' num2str(psver) '.mat'];
apsbandssbname = ['./tca_bands_sb' num2str(psver) '.mat'];
apsname=['./tca',num2str(psver)];
apssbname=['./tca_sb',num2str(psver)];
iononame = ['./ica' num2str(psver) '.mat'];
ionosbname = ['./ica_sb' num2str(psver) '.mat'];
tidename=['./tide',num2str(psver)];
tidesbname=['./tide_sb',num2str(psver)];
hgtname=['./hgt',num2str(psver)];


sclaname=['./scla',num2str(psver)];
sclasbname=['./scla_sb',num2str(psver)];
sclasmoothname=['./scla_smooth',num2str(psver)];
sclasbsmoothname=['./scla_smooth_sb',num2str(psver)];
meanvname=['./mv',num2str(psver)];
meanVname=['./mean_v.mat'];
ztdname=['./gacos.mat'];

drop_ifg_index=getparm('drop_ifg_index');
small_baseline_flag=getparm('small_baseline_flag');
scla_deramp = getparm('scla_deramp');
fig_name_suffix = '';
bands = [];
unwrap_ifg_index=[1:ps.n_image];


if unwrap_ifg_index(1)~=ps.master_ix & unwrap_ifg_index(end)~=ps.master_ix
    unwrap_ifg_index=setdiff(unwrap_ifg_index,ps.master_ix); % need to include it if not ifgs either side of master
end
if ~isempty(ifg_list)
     unwrap_ifg_index=intersect(unwrap_ifg_index,ifg_list);
     ifg_list=[];
end


if strcmpi(small_baseline_flag,'y')
    phuwres=load(phuwsbresname,'sm_cov');
    if isfield(phuwres,'sm_cov');
         sm_cov=phuwres.sm_cov(unwrap_ifg_index,unwrap_ifg_index);
    else
         sm_cov=eye(length(unwrap_ifg_index));
    end
else
    
    if ~exist([ifgstdname,'.mat',],'file')
    sm_cov=eye(length(unwrap_ifg_index));
    else
        ifgstd=load(ifgstdname);
        if isfield(ifgstd,'ifg_std');
            ifgvar=(ifgstd.ifg_std*pi/181).^2;
            sm_cov=diag(ifgvar(unwrap_ifg_index));
        else
            sm_cov=eye(length(unwrap_ifg_index));
        end
    end
end


day=day(unwrap_ifg_index,:);
G=[ones(size(day)),day-master_day]; 
lambda=getparm('lambda');
ph_all=zeros(n_ps,1);
textsize=0;
units='mm/yr';

if isempty(ifg_list)
    ifg_list=1:size(ph_all,2);
end
n_ifg_plot=length(ifg_list);


group_type=value_type;

if strcmpi(group_type,'meanv')
      disp('load data ...')
      uw=load(phuwname);
      ph_uw=double(uw.ph_uw);
      in = inpolygon(lonlat(:,1),lonlat(:,2),plgn(:,1),plgn(:,2));
      scla=load(sclaname);
      if isfield(scla,'C_ps_uw')
          ph_uw=ph_uw - scla.ph_scla - repmat(scla.C_ps_uw,1,size(ph_uw,2)); % master phase doesn't effect plot, but better for ts plot
      else
          ph_uw=ph_uw - scla.ph_scla;
      end
      ph_uw     =ph_uw(:,unwrap_ifg_index);
      ref_ps    =ps_setref;
      ph_uw     =ph_uw-repmat(nanmean(ph_uw(ref_ps,:),1),n_ps,1);
      clear uw scla
      fig_name  = 'detrend_meanV';
      m         = lscov(G,double(ph_uw'),sm_cov);
      meanV     = -m(2,:)'*365.25/4/pi*lambda*1000; 
      lonlat_in = lonlat(in,:);
      meanV_in  = meanV(in,:);
      fit_pl    = fit([lonlat_in(:,1), lonlat_in(:,2)], meanV_in,'poly11');
      trend     = fit_pl.p10 * lonlat(:,1) +  fit_pl.p01 * lonlat(:,2) + fit_pl.p00;
      ph_disp   = meanV - trend;
      save ps_plot_dt ph_disp  meanV trend plgn '-v7.3' 
      
  elseif strcmpi(group_type, 'rmeq')
      disp('load data ...')
      uw=load(EQname);
      ph_uw=double(uw.ph_uw);
      in = inpolygon(lonlat(:,1),lonlat(:,2),plgn(:,1),plgn(:,2));
      ph_uw=ph_uw(:,unwrap_ifg_index);
      ref_ps=ps_setref;
      ph_uw=ph_uw-repmat(nanmean(ph_uw(ref_ps,:),1),n_ps,1);
      clear uw scla
      fig_name = 'detrend_meanV_rmEQ';
      m=lscov(G,double(ph_uw'),sm_cov);
      ph_all=-m(2,:)'*365.25/4/pi*lambda*1000; 
      lonlat_in = lonlat(in,:);
      ph_all_in = ph_all(in,:);
      fit_pl = fit([lonlat_in(:,1), lonlat_in(:,2)],ph_all_in,'poly11');
      model = fit_pl.p10 * lonlat(:,1) +  fit_pl.p01 * lonlat(:,2) + fit_pl.p00;
      ph_all=ph_all-model;
      ph_disp=ph_all(:,ifg_list);
      save ps_plot_dtrnd_rmEQ ph_disp '-v7.3'

  elseif strcmpi(group_type, 'gps')
      fig_name = 'detrend_meanV';
      disp('load data ...')
      uw=load(phuwname);
      ph_uw=double(uw.ph_uw);
      scla=load(sclaname);
      if isfield(scla,'C_ps_uw')
          ph_uw=ph_uw - scla.ph_scla - repmat(scla.C_ps_uw,1,size(ph_uw,2)); % master phase doesn't effect plot, but better for ts plot
      else
          ph_uw=ph_uw - scla.ph_scla;
      end
      clear uw scla
      ph_uw     = ph_uw(:,unwrap_ifg_index);
      m         = lscov(G,double(ph_uw'),sm_cov);
      ph_all    = -m(2,:)'*365.25/4/pi*lambda*1000; 
      inPS  = inpolygon(lonlat(:,1),lonlat(:,2),plgn(:,1),plgn(:,2));
      inGPS = inpolygon(Gdat(:,1),Gdat(:,2),plgn(:,1),plgn(:,2));
      if hdng > -45
        LosGps    = gps2los([Gdat(:,1) Gdat(:,2) Gdat(:,3)*-1 Gdat(:,4)]);
        else
        LosGps    = gps2los(Gdat(inGPS,1:4));
      end      
      
      if ~isempty(ph_all(inPS))
         ph_all = ph_all - (nanmean(ph_all(inPS)) - nanmean(LosGps(:,3)));
      else
         ph_all = ph_all - nanmean(LosGps(:,3));
      end
      ph_disp=ph_all(:,ifg_list);
      save ps_plot_dtrnd ph_disp '-v7.3'  
      
  elseif strcmpi(group_type, 'inShpGps')
      disp('load data .......')
      fig_name  = 'detrend_meanV';
      uw        = load(phuwname);
      ph_uw     = double(uw.ph_uw);
      scla      = load(sclaname);
      if isfield(scla,'C_ps_uw')
          ph_uw = ph_uw - scla.ph_scla - repmat(scla.C_ps_uw,1,size(ph_uw,2)); % master phase doesn't effect plot, but better for ts plot
      else
          ph_uw = ph_uw - scla.ph_scla;
      end
      clear uw scla
      %%% calculate mean Velocity
      ph_uw     = ph_uw(:,unwrap_ifg_index);
      m         = lscov(G,double(ph_uw'),sm_cov);
      meanV     = -m(2,:)'*365.25/4/pi*lambda*1000; 
      %%% model GPS
      inGPS2    = inpolygon(Gdat(:,1),Gdat(:,2),plgn(:,1),plgn(:,2));
      GPS4Model = Gdat(inGPS2,:);
      modelGPS  = gps2los(lonlat,GPS4Model);
      rsdual    = meanV - modelGPS;
      %%% extract trend
      fit_pl    = fit([lonlat(:,1), lonlat(:,2)], rsdual, 'poly11');
      trend     = feval(fit_pl,[lonlat(:,1), lonlat(:,2)]);
      %%% save  
      ph_disp   = meanV(:,ifg_list) - trend;
      save ps_plot_dtrnd ph_disp '-v7.3'
      save gps2losModel lonlat meanV trend modelGPS GPS4Model '-v7.3'

  elseif strcmpi(group_type, 'inFrmGps')
      disp('load data .......')
      fig_name  = 'detrend_meanV';
      uw        = load(phuwname);
      ph_uw     = double(uw.ph_uw);
      scla      = load(sclaname);
      if isfield(scla,'C_ps_uw')
          ph_uw = ph_uw - scla.ph_scla - repmat(scla.C_ps_uw,1,size(ph_uw,2)); % master phase doesn't effect plot, but better for ts plot
      else
          ph_uw = ph_uw - scla.ph_scla;
      end
      clear uw scla
      %%% calculate mean Velocity
      ph_uw     = ph_uw(:,unwrap_ifg_index);
      m         = lscov(G,double(ph_uw'),sm_cov);
      meanV     = -m(2,:)'*365.25/4/pi*lambda*1000; 
      %%% model GPS
      GPS4Model = Gdat(inFrmGPS,:);
      modelGPS  = gps2los(lonlat,GPS4Model);
      rsdual    = meanV - modelGPS;
      %%% extract trend
      fit_pl    = fit([lonlat(:,1), lonlat(:,2)], rsdual, 'poly11');
      trend     = feval(fit_pl,[lonlat(:,1), lonlat(:,2)]);
      %%% save  
      ph_disp   = meanV(:,ifg_list) - trend;
      save ps_plot_dtrnd ph_disp '-v7.3'
      save gps2losModel lonlat meanV trend modelGPS GPS4Model '-v7.3'
      
  elseif strcmpi(group_type, 'IntSeis')
      disp('load data .......')
      fig_name = 'detrend_meanV';
      uw       = load(phuwname);
      ph_uw    = double(uw.ph_uw);
      scla     = load(sclaname);
      if isfield(scla,'C_ps_uw')
         ph_uw = ph_uw - scla.ph_scla - repmat(scla.C_ps_uw,1,size(ph_uw,2)); % master phase doesn't effect plot, but better for ts plot
      else
         ph_uw = ph_uw - scla.ph_scla;
      end
      clear uw scla
      %%% calculate mean Velocity
      ph_uw    = ph_uw(:,unwrap_ifg_index);
      m        = lscov(G,double(ph_uw'),sm_cov);
      meanV    = -m(2,:)'*365.25/4/pi*lambda*1000; 
      %%% load interseismic model
      ISmdl    = load('');
      inFrmIS  = inpolygon(ISmdl(:,1), ISmdl(:,2), frm(:,1), frm(:,2));
      modelIS  = gps2los(lonlat,inFrmIS(:,1:4));  % model must be in this order; Longitude, Latitude, East, North
      rsdl     = meanV - modelIS;
      %%% extract trend
      fit_pl   = fit([lonlat(:,1), lonlat(:,2)], rsdl, 'poly11');
      trend    = feval(fit_pl,[lonlat(:,1), lonlat(:,2)]);
      %%% save  
      ph_disp   = meanV(:,ifg_list) - trend;
      save ps_plot_dtrnd ph_disp '-v7.3'
      save IntSeisModel lonlat meanV trend modelIS '-v7.3'
      
  elseif strcmpi(group_type, 'gpslos2')
      %%% get GPS
      disp('load GPS .......')
      fltM = ['/nas3/tohid/stamps_result/malazgirt/others/fault_modl.mat'];
      if exist(fltM, 'file')
      flt=load(fltM);
      fLon=flt.flt(:,1); fLat=flt.flt(:,2);
      indF = find(fLon>min(lonlat(:,1)) & fLon<max(lonlat(:,1)));
      plot(fLon(indF),fLat(indF),'k','linewidth',1.5)
      end
      [pLon, pLat] = getpts;
      frm = [pLon  pLat; pLon(1) pLat(1)];
      fig_name = 'detrend_meanV';
      Gdat = load('gps.dat');
      plot(frm(:,1),frm(:,2),'r')
      inGPS = inpolygon(Gdat(:,1),Gdat(:,2),frm(:,1),frm(:,2));
      quiver(Gdat(inGPS,1),Gdat(inGPS,2),Gdat(inGPS,3),Gdat(inGPS,4))
      plot(Gdat(inGPS,1),Gdat(inGPS,2),'.r','markersize',10)
      pause(3)
      %%% load data
      disp('load data .......')
      uw=load(phuwname);
      ph_uw=double(uw.ph_uw);
      scla=load(sclaname);
      if isfield(scla,'C_ps_uw')
          ph_uw=ph_uw - scla.ph_scla - repmat(scla.C_ps_uw,1,size(ph_uw,2)); % master phase doesn't effect plot, but better for ts plot
      else
          ph_uw=ph_uw - scla.ph_scla;
      end
      clear uw scla
      ph_uw    = ph_uw(:,unwrap_ifg_index);
      m        = lscov(G,double(ph_uw'),sm_cov);
      ph_all   = -m(2,:)'*365.25/4/pi*lambda*1000; 
      modelGPS = gps2los(lonlat,Gdat(inGPS,:));
      ph_all   = ph_all - modelGPS;
      %%% remove trend
      in = inpolygon(lonlat(:,1),lonlat(:,2),plgn(:,1),plgn(:,2));
      if length(lonlat(in,:)) ~= 0
        lonlat_in = lonlat(in,:); ph_all_in = ph_all(in,:);
        fit_pl = fit([lonlat_in(:,1), lonlat_in(:,2)],ph_all_in,'poly11');
        model = feval(fit_pl,[lonlat(:,1), lonlat(:,2)]);
        ph_all=ph_all-model;
      end
      ph_disp=ph_all(:,ifg_list);
      save ps_plot_dtrnd ph_disp '-v7.3' 
      save nearfltgps lonlat ph_disp  modelGPS '-v7.3'
      
  elseif strcmpi(group_type, 'losgps')
      fig_name = 'detrend_meanV';
      disp('load data ...')
      uw=load(phuwname);
      ph_uw=double(uw.ph_uw);
      scla=load(sclaname);
      if isfield(scla,'C_ps_uw')
          ph_uw=ph_uw - scla.ph_scla - repmat(scla.C_ps_uw,1,size(ph_uw,2)); % master phase doesn't effect plot, but better for ts plot
      else
          ph_uw=ph_uw - scla.ph_scla;
      end
      clear uw scla
      ph_uw     = ph_uw(:,unwrap_ifg_index);
      ref_ps    = ps_setref;
      ph_uw     = ph_uw-repmat(nanmean(ph_uw(ref_ps,:),1),n_ps,1);
      m         = lscov(G,double(ph_uw'),sm_cov);
      ph_los  = -m(2,:)'*365.25/4/pi*lambda*1000; 
      
      %%% convert los to E-W direction
      hdng    = getparm('heading');
      load fault_strike
      strk = extFaultAz;
      la   = extLA;
      ph_all = ph_los ./ (sind(la).*(cosd(hdng).*sind(strk)-sind(hdng).*cosd(strk)));
      %%% fit to the E-W velocity field of GPS
      if hdng > -45
        fit_GPS   = fit([Gdat(:,1), Gdat(:,2)],Gdat(:,3)*-1,'poly11');
        else
        fit_GPS   = fit([Gdat(:,1), Gdat(:,2)],Gdat(:,3),'poly11');
      end     
      modelGPS  = fit_GPS.p10 * lonlat(:,1) +  fit_GPS.p01 * lonlat(:,2) + fit_GPS.p00;
      ph_all_in = ph_all - modelGPS;
      fit_PS    = fit([lonlat(:,1), lonlat(:,2)],ph_all,'poly11');
      modelPS   = fit_PS.p10 * lonlat(:,1) +  fit_PS.p01 * lonlat(:,2) + fit_PS.p00;
      ph_all    = ph_all - modelPS;
      %%% reference to polygon 
      inPS  = inpolygon(lonlat(:,1),lonlat(:,2),plgn(:,1),plgn(:,2));
      fit_PS    = fit([lonlat(inPS,1), lonlat(inPS,2)],ph_all(inPS),'poly11');
      modelPS   = fit_PS.p10 * lonlat(:,1) +  fit_PS.p01 * lonlat(:,2) + fit_PS.p00;
      ph_all    = ph_all - modelPS;
      ph_disp=ph_all(:,ifg_list);
      save ps_plot_dtrnd ph_disp '-v7.3'
      
  elseif strcmpi(group_type, 'phase')
      uw=load(phuwname);
      ph_uw=double(uw.ph_uw);
      ph_uw=ph_uw(:,unwrap_ifg_index);
      clear uw
      in = inpolygon(lonlat(:,1),lonlat(:,2),plgn(:,1),plgn(:,2));
      lonlat_in = lonlat(in==1,:);
      ph_uw_in = ph_uw(in==1,:);
      dtrnd_phuw = [];
        for i=1:size(ph_uw,2)
           fit_pl = fit([lonlat_in(:,1), lonlat_in(:,2)],ph_uw_in(:,i),'poly11');
           model = fit_pl.p10 * lonlat(:,1) +  fit_pl.p01 * lonlat(:,2) + fit_pl.p00;
           dtrnd_phuw(:,i) = ph_uw(:,i)-model;
        end
      scla=load(sclaname);
      if isfield(scla,'C_ps_uw')
          ph_uw=ph_uw - scla.ph_scla(:,unwrap_ifg_index) - repmat(scla.C_ps_uw,1,size(ph_uw,2)); % master phase doesn't effect plot, but better for ts plot
      else
          ph_uw=ph_uw - scla.ph_scla(:,unwrap_ifg_index);
      end
      save phuw_dtrnd ph_uw '-v7.3'
      clear uw scla
      fig_name = 'detrend_phase';
      G=[ones(size(day)),day-master_day]; 
      lambda=getparm('lambda');
      m=lscov(G,double(ph_uw'),sm_cov);
      ph_all=-m(2,:)'*365.25/4/pi*lambda*1000; % m(1,:) is master APS + mean deviation from model
      ph_disp=ph_all(:,ifg_list);
  else
            error('unknown value type')
end












xgap=0.1;
ygap=0.2;
[Y,X]=meshgrid([0.7:-1*ygap:0.1],[0.1:xgap:0.8]);

if ~isempty(lon_rg)
    ix=lonlat(:,1)>=lon_rg(1)&lonlat(:,1)<=lon_rg(2);
    lonlat=lonlat(ix,:);
end

if ~isempty(lat_rg)
    ix=lonlat(:,2)>=lat_rg(1)&lonlat(:,2)<=lat_rg(2);
    lonlat=lonlat(ix,:);
end

max_xy=llh2local([max(lonlat),0]',[min(lonlat),0]);
fig_ar=4/3; % aspect ratio of figure window
useratio=1; % max fraction of figure window to use
n_i=max_xy(2)*1000;
n_j=max_xy(1)*1000;
ar=max_xy(1)/max_xy(2); % aspect ratio (x/y)
if n_x==0
    n_y=ceil(sqrt((n_ifg_plot)*ar/fig_ar)); % number of plots in y direction
    n_x=ceil((n_ifg_plot)/n_y);
    fixed_fig = 1; 			% figure with fixed aspect ratio
else
    n_y=ceil((n_ifg_plot)/n_x);
    fixed_fig = 0;
end


d_x=useratio/n_x;
d_y=d_x/ar*fig_ar;
if d_y>useratio/n_y & fixed_fig==1	% TS figure exceeds fig size
    d_y=useratio/n_y; 
    d_x=d_y*ar/fig_ar;
    h_y=0.95*d_y;
    h_x=h_y*ar/fig_ar;

    fig_size=0;
elseif d_y>useratio/n_y & fixed_fig==0 
    h_y=0.95*d_y;
    h_x=h_y*ar/fig_ar;

    y_scale = d_y*n_y;
    d_y=d_y/y_scale;   
    fig_size=1;		% check to indicate fig needs to be adapted
    h_y=0.95*d_y;

else
    h_y=0.95*d_y;
    h_x=h_y*ar/fig_ar;
    fig_size=0;
end
y=1-d_y:-d_y:0;
x=1-useratio:d_x:1-d_x;


[imY,imX]=meshgrid(y,x);
if textsize==0
    textsize=round(10*4/n_x);
    if textsize>16
        textsize=16;
    elseif textsize<8
        textsize=8;
    end
end

l_t=1/9*abs(textsize)/10; % text length
h_t=1/50*abs(textsize)/10; % text height
x_t=round((h_x-l_t)/h_x/2*n_j);
y_t=round(h_t*1.2/h_y*n_i);


if isreal(ph_all)
    if ref_ifg~=0
        if ref_ifg==-1
            ph_disp=ph_disp-[ph_disp(:,1),ph_disp(:,1:end-1)];
        else
            ph_disp=ph_disp-repmat(ph_all(:,ref_ifg),1,size(ph_disp,2));
        end
    else
        ref_ifg=master_ix;
    end
    if ref_ps~=0
        ref_ph=(ph_disp(ref_ps,:));
        mean_ph=zeros(1,size(ph_disp,2));
        for i=1:size(ph_disp,2)
            mean_ph(i)=mean(ref_ph(~isnan(ref_ph(:,i)),i),1);
            if isnan(mean_ph(i))
                mean_ph(i)=0;
                fprintf(['Interferogram (' num2str(i) ') does not have a reference area\n'])
            end
        end
        clear i
        ph_disp=ph_disp-repmat(mean_ph,n_ps,1);
    end
    phsort=sort(ph_disp(~isnan(ph_disp)));
    
    if isempty(lims)
        maxph=phsort(round(length(phsort)*.999));
        minph=phsort(ceil(length(phsort)*.001));
        lims=[minph,maxph];
        if isempty(lims)==1
            fprintf(['Interferograms do not contain data.\n'])
        end
        
    end
else
    if ref_ifg==0
        ref_ifg=master_ix;
    elseif ref_ifg==-1
        ph_disp=ph_disp.*conj([ph_disp(:,1),ph_disp(:,1:end-1)]);
    end
    if ref_ps~=0
        ph_disp=ph_disp./abs(ph_disp);
        ref_ph=(ph_disp(ref_ps,:));
        mean_ph=zeros(1,size(ph_disp,2));
        for i=1:size(ph_disp,2)
            mean_ph(i)=sum(ref_ph(~isnan(ref_ph(:,i)),i));
        end
        clear i
        ph_disp=ph_disp.*conj(repmat(mean_ph,n_ps,1));
    end
    lims=[-pi,pi];

end

if plot_flag==-1
    %savename=['~/ps_plot_',value_type]
    h_fig=[];
    lims=[];
    ifg_data_RMSE=[];
    h_axes_all=[];
    savename=['ps_plot_',value_type]
    try
       save(savename,'ph_disp','ifg_list','-v7.3')
    catch
       save(['~/',savename],'ph_disp','ifg_list','-v7.3')
       fprintf('Warning: Read access only, values in home directory instead\n')
    end
else
  h_fig = figure;
  set(gcf,'renderer','zbuffer','name',fig_name)

  if fig_size==1
      Position = get(gcf,'Position');
      Position(1,2)=50;
      Position(1,4)=Position(1,2)+Position(1,4)*y_scale;
      set(gcf,'Position',Position)
  end


  i_im=0;
  ifg_data_RMSE=NaN([size(ph_disp,2) 1]) ;
  if size(ifg_list,1)>1
    ifg_list = ifg_list';
  end
  for i=ifg_list
    i_im=i_im+1;
    if n_ifg_plot>1
        h_axes = axes('position',[imX(i_im),imY(i_im),h_x,h_y]);
        set(h_axes,'Xcolor',[1 1 1]);			% [DB] remove the contours of the axes figure
        set(h_axes,'Ycolor',[1 1 1]);
        
        h_axes_all(i)=h_axes;
        clear h_axes
    else
        h_axes_all(i)=gca;
    end
    % check if external data is requested to be plotted
    if ext_data_flag==1
        if exist(ext_data_path,'file')==2
            % this is an individual file identical for each ifgs
            loadname=ext_data_path;
        else      
            if size(day,1)==size(ph_all,2) && strcmpi(small_baseline_flag,'y')
                % This is for external data for a SB to a SM inverted network
                if aps_band_flag==1
                    fprintf('Still to do')
                    keyboard 
                else                   
                    % checking if there is data for this interferogram
                    loadname = [ext_data_path filesep  datestr(ps.day(i,1),'yyyymmdd') '_' datestr(ps.master_day,'yyyymmdd') '.mat'];
                end
 
            else
                % This is external data for a SM or SB network. It is not a
                % SB to SM inverted network.
                % The path is specified
                if aps_band_flag==1
                    % checking if there is data for this interferogram
                    loadname = [ext_data_path filesep datestr(ps.ifgday(ifg_number,1),'yyyymmdd') '_' datestr(ps.ifgday(ifg_number,2),'yyyymmdd') '.mat'];
                else
                    % checking if there is data for this interferogram
                    loadname = [ext_data_path filesep datestr(ps.ifgday(i,1),'yyyymmdd') '_' datestr(ps.ifgday(i,2),'yyyymmdd') '.mat'];
                end
 
            end
            
            
            
           
        end
        if exist(loadname,'file')==2
             % loading of the data 
             ext_data = load(loadname);
             
             
             % Checking if there is a ph_disp variable stored in the files 
             if isfield(ext_data,{'ph_disp'})
                 % Correct the data with respect to the minimization of the
                 % residual between this dataset and the data
                 % convert to a local reference
                 n_data_points = size(ext_data.lonlat,1);
                 data_xy=llh2local(ext_data.lonlat',ps.ll0)*1000;
                 ps_xy=llh2local(ps.lonlat',ps.ll0)*1000;

                 dist_sq=(repmat(ps_xy(1,:)',1,n_data_points) -repmat(data_xy(1,:),ps.n_ps,1)).^2+(repmat(ps_xy(2,:)',1,n_data_points) -repmat(data_xy(2,:),ps.n_ps,1)).^2; 
                 for kk=1:n_data_points      
                     ref_data= (dist_sq(:,kk)<=ref_radius_data^2);
                     if sum(ref_data)>0
                        ifg_mean_data_point(kk,1) = mean(ph_disp(ref_data,i_im));
                     else
                        ifg_mean_data_point(kk,1) = NaN;  
                     end
                 end
                 clear ref_data
                 if ~isempty(ifg_mean_data_point) && sum(isnan(ifg_mean_data_point))~=n_data_points
                     % mean residual 
                     ext_data.mean_residual = mean(ext_data.ph_disp(~isnan(ifg_mean_data_point),1)-ifg_mean_data_point(~isnan(ifg_mean_data_point)));
                     % correct data phases
                     ext_data.ph_disp(:,1) = ext_data.ph_disp(:,1)-ext_data.mean_residual;
                     % computing the RMSE
                     ext_data.residual = ext_data.ph_disp(~isnan(ifg_mean_data_point),1)-ifg_mean_data_point(~isnan(ifg_mean_data_point));
                     ext_data.RMSE = sqrt(mean((ext_data.residual).^2));
                     clear ifg_mean_data_point
                     ifg_data_RMSE(i)=ext_data.RMSE;
                 else
                    fprintf(['No observation within the ref_radius_data, increase the size. \n']) 
                    fprintf(['External data not plotted for this interferogram. \n']) 
                    ext_data = [];
                 end
                 clear mean_data_point_residual

    %              % correct the data with respect to the reference area
    %              point_ref = ps_setref(ext_data);
    %              ext_data.ph_disp(:,1) = ext_data.ph_disp(:,1)-mean(ext_data.ph_disp(point_ref,1));
            end

        else
           ext_data = []; 
        end
    else
       ext_data = []; 
    end
    
    
    ps_plot_ifg(ph_disp(:,i_im),plot_flag,lims,lon_rg,lat_rg,ext_data);
    %plot_phase(ph_tc(:,i)*conj(ph_tc(ref_ix,i)));
    box on
    if n_ifg_plot>1
        set(gca,'yticklabel',[])
        set(gca,'xticklabel',[])
    end
    xlim=get(gca,'xlim');
    x_t=(h_x-l_t)/2/h_x*(xlim(2)-xlim(1))+xlim(1);
    ylim=get(gca,'ylim');
    if textsize>0
        y_t=(h_y-1.2*h_t)/h_y*(ylim(2)-ylim(1))+ylim(1);
    else
        y_t=(0.5*h_t)/h_y*(ylim(2)-ylim(1))+ylim(1);
    end
    %xlabel([num2str((day(i)/365.25),3),'yr, ',num2str(round(bperp(i))),'m'])
    
        if textsize~=0 & size(day,1)==size(ph_all,2) & aps_band_flag==0 && strcmpi(small_baseline_flag,'n')
            % text for SM ifgs
            t=text(x_t,y_t,[datestr(day(i),'dd mmm yyyy')]);
            set(t,'fontweight','bold','color',textcolor,'fontsize',abs(textsize))
        elseif textsize~=0 & ps.n_ifg==size(ph_all,2) & aps_band_flag==0 && strcmpi(small_baseline_flag,'y')
            % text for SB ifgs
            t=text(x_t,y_t,['      ifg ' num2str(i)]);
            set(t,'fontweight','bold','color',textcolor,'fontsize',abs(textsize))
        elseif textsize~=0 & size(bands,1)==size(ph_all,2) & aps_band_flag==1
            % text for band filtered data
            t=text(x_t,y_t,[num2str(round((bands(i,1)/100))*100/1000) ' - ' num2str(round((bands(i,2)/100))*100/1000) ' km']);
            set(t,'fontweight','bold','color',textcolor,'fontsize',abs(textsize))
        end
    
    
    if cbar_flag==0 & (i==ref_ifg | (isempty(intersect(ref_ifg,ifg_list)) & i==ifg_list(1))) 
        if n_ifg_plot>1
            
           
            h=colorbar('South');
            xlim=get(h,'xlim');
            set(h,'xlim',[xlim(2)-64,xlim(2)])

        else
            %h=colorbar('SouthOutside');
            h = colorbar('peer',gca);
            ylim=get(h,'ylim');
            set(h,'ylim',[ylim(2)-64,ylim(2)])
        end

        if diff(lims)>1 | diff(lims)==0
            plotlims=round(lims*10)/10;
        else
            limorder=ceil(-log10(diff(lims)))+2;
            plotlims=round(lims*10^limorder)/10^limorder;
        end
        if n_ifg_plot>1
                set(h,'xtick',[xlim(2)-64,xlim(2)],'Xticklabel',plotlims,'xcolor',textcolor,'ycolor',textcolor,'fontweight','bold','color',textcolor,'FontSize',abs(textsize))
                h=xlabel(h,units);
                pos=get(h,'position');
                pos(2)=pos(2)/2.2;
                set(h,'position',pos,'FontSize',abs(textsize));
        else
                set(h,'ytick',[ylim(2)-64,ylim(2)],'yticklabel',plotlims,'xcolor',textcolor,'ycolor',textcolor,'fontweight','bold','color',textcolor,'FontSize',abs(textsize))
                set(get(h,'ylabel'),'String',units,'FontSize',abs(textsize))  

        end
    end
  end
end

if ext_data_flag==1
    for k=1:size(ifg_data_RMSE,1)
        if k==1
            fprintf(['\nRMSE between interferogram(s) and external data \n'])
        end
        
        if aps_band_flag==1
            fprintf(['Band' num2str(k) ' : ' num2str(ifg_data_RMSE(k)) '\n']);
        else
            if size(ifg_data_RMSE,1)==ps.n_ifg
                fprintf(['ifg ' num2str(k) ' \t ' datestr(ps.ifgday(k,1),'yyyymmdd') '-' datestr(ps.ifgday(k,2),'yyyymmdd') ' \t ' num2str(ifg_data_RMSE(k)) '\n']);
            else
                 fprintf(['dataset ' num2str(k) ' \t ' num2str(ifg_data_RMSE(k)) '\n']);               
            end
        end
    end
end
hold on
if exist('point.dat','file') == 2
 fid=fopen('point.dat');
  N2=textscan(fid,'%f %f',[2 inf]);
       x4=N2{1};
       y4=N2{2};
       scatter(x4,y4,'filled','black')
       %text(x4,y4,n4)
       hold on
end
if exist('plot_color_scheme_old','var')==1
   setparm('plot_color_scheme',plot_color_scheme_old) 
end

fprintf('Color Range: %g to %g %s\n',lims,units)
if ts_flag == 1
  figure(h_fig);
  clear all % clean up to save memory
  
  % Place button for new TS plots
  fPosition=get(gcf,'position');
  %  pos_new = [pos(1) pos(2)+0.1 pos(3) pos(4)-0.1]
  
  % new TS plot button
  mButton=uicontrol('Style', 'pushbutton', 'Callback', 'clear ph_uw; ts_plot',...
    'String','TS plot', 'Position', [150 30 90 20] ); % if exist don't create
  %                                      left bottom width height
  %mButtonposition=get(mButton,'Position')
 
  % new TS plot button for selecting two points and difference
  mButton=uicontrol('Style', 'pushbutton', 'Callback', 'clear ph_uw; ts_plotdiff',...
    'String','TS double diff.', 'Position', [470 30 90 20] ); % 
  
  % Radius Text boxes
  mTextBox=uicontrol('Style', 'edit','String','radius (m) ', 'Position',...
      [320 30 90 20] );

  mEditBox=uicontrol('Style', 'edit','String','100', 'Position',...
      [410 30 30 20],'BackgroundColor',[1 1 1] );
%  ts_plot   % select a point then plot, for the first time.
end

save('prof_dtrnd','lonlat','ph_disp','-v7.3')
%%% plot mean velocity profile
pushbutton=uicontrol('Style', 'pushbutton', 'Callback', 'tps_profile_dtrnd',...
    'String','Profile', 'Position', [20 50 75 20] ); % 
%%% plot displacement profile
pushbutton=uicontrol('Style', 'pushbutton', 'Callback', 'dpl_profile_dtrnd',...
    'String','Displacement', 'Position', [20 75 75 20] ); % 
    
TextBox=uicontrol('Style', 'edit','String','100', 'Position',...
      [20 25 50 20],'BackgroundColor',[1 1 1] );


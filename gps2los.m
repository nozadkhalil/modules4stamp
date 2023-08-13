function model = gps2los(lonlat,GPS)

%%%  GPS file input format; [lon lat E-W N-S Up]
%%%  lonlat is [lon lat] of los 
%%%  this script fit a second degree surface to E-W and N-S velocities of GPS 
%%%  then extract the value of fit for each pixel and convert them to los


disp('calculate model .......')
%%% LOAD UNIT VECTOR
%%% order of unite vectors [lon lat elv E-W N-S Up]
UV=load('unit_vector.ll');
UV=UV(1:2:end,:);

%%% FIT SURFACE to unit vectors 2nd degree for X and Y axis
UVE  = fit([UV(:,1), UV(:,2)], UV(:,4),'poly22');
UVN  = fit([UV(:,1), UV(:,2)], UV(:,5),'poly22');

%%% caculate unit vector value for each pixel
ell  = feval(UVE,[lonlat(:,1), lonlat(:,2)]);
nll  = feval(UVN,[lonlat(:,1), lonlat(:,2)]);

%%% FIT SURFACE to GPS 2nd degree for X and Y axis
fite  = fit([GPS(:,1), GPS(:,2)],GPS(:,3),'poly22');
fitn  = fit([GPS(:,1), GPS(:,2)],GPS(:,4),'poly22');

%%% caculate model value for each pixel
mode  = feval(fite,[lonlat(:,1), lonlat(:,2)]);
modn  = feval(fitn,[lonlat(:,1), lonlat(:,2)]);

%%% CONVERT GPS to line of sight direction
model = (ell.*mode)+(nll.*modn); 



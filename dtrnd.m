function dtrnd
disp('Remove the Trend')
TextBox2=findobj('Background',[1 1 1]);
value_type = get(TextBox2,'String');					
value_type = lower(value_type);

if strcmpi(value_type,'meanv')
   dtrnd_psplot('meanV')
elseif strcmpi(value_type,'rmeq')
   dtrnd_psplot('rmeq')
elseif strcmpi(value_type,'gps')
   dtrnd_psplot('gps')
elseif strcmpi(value_type,'losgps')
   dtrnd_psplot('losgps')
elseif strcmpi(value_type,'gpslos2')
   dtrnd_psplot('gpslos2')
elseif strcmpi(value_type,'inShpGps')
   dtrnd_psplot('inShpGps')
elseif strcmpi(value_type,'inFrmGps')
   dtrnd_psplot('inFrmGps')
else
   dtrnd_psplot('phase')
end

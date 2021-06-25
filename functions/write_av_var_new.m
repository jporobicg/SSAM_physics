% WRITE_BOX_AV
%
%  tims    Base [1990 1 1 +10], so tim = tims - greg2time([1990 1 1 10 0 0]);
%  bxid    If not just 0:(nbox-1)
%  fnm     output file name, incl '.nc'.  [default 'av_prop.nc']
%  vars    vector (ascending order):  1=temp  2=salt  3=w_flux  4=bx_nett
%                         5=eta
%  dat     cell array of matrices corresponding to 'vars'
%
% Jeff Dunn CMAR  19/1/06
%
% USAGE: write_box_av(tims,bxid,fnm,vars,dat)
%tims=nctime;
%temp=temperature;
%salt=salinity;
%w= vertical;


function write_av_var_new(tims, bid, temp, salt, w, fnm)

nbx = length(bid);
ntm = length(tims);
nlv = size(temp,3);


% if min(tims) > 10000
%    t90  = greg2time([1990 1 1 10 0 0]);
%    tims = tims-t90;
% end

%nc = netcdf(fnm,'noclobber');
nc_new = netcdf.create(fnm, '64BIT_OFFSET');

% --- Create global attributes


% dimension variables
time      =  netcdf.defDim(nc_new, 'time', netcdf.getConstant('NC_UNLIMITED'));
level     =  netcdf.defDim(nc_new, 'level', nlv);
boxes     =  netcdf.defDim(nc_new, 'boxes', nbx);
NC_GLOBAL = netcdf.getConstant('NC_GLOBAL'); % constant
%% Variables and attributes
var_w          =  netcdf.defVar(nc_new, 'verticalflux', 'float', [level  boxes time]);
var_salinity   =  netcdf.defVar(nc_new, 'salinity', 'float', [level boxes time]);
var_temp       =  netcdf.defVar(nc_new, 'temperature', 'float', [level boxes time]);
time           =  netcdf.defVar(nc_new, 'time', 'double', time);
boxes           =  netcdf.defVar(nc_new, 'boxes', 'int', boxes);
level           =  netcdf.defVar(nc_new, 'level', 'int', level);
%ncdim('time', 0, nc);
%ncdim('level', nlv, nc);
%ncdim('boxes', nbx, nc);
% time in seconds
%% time_from
netcdf.putAtt(nc_new, time , 'long_name', char('simulation_time'));
netcdf.putAtt(nc_new, time , 'units', char('seconds since 1970-01-01 00:00:00'));
netcdf.putAtt(nc_new, time , 'grads_step', char('seconds'));
netcdf.putAtt(nc_new, time , 'dt', 21600);
%% Boxes
netcdf.putAtt(nc_new, boxes , 'long_name', char('Box IDs'));
%% Levels
netcdf.putAtt(nc_new, level , 'long_name', char('layer index; 1=near-surface'));
%% Vertical flux
%% W - velocity component
netcdf.putAtt(nc_new, var_w , 'long_name', char('ocean vertical velocity'));
netcdf.putAtt(nc_new, var_w , 'units', char('meter second-1'));
netcdf.putAtt(nc_new, var_w , 'field', char('v-velocity,  scalar, series'));
%% Salinity
netcdf.putAtt(nc_new, var_salinity , 'long_name', char('sea_water_reference_salinity'));
netcdf.putAtt(nc_new, var_salinity , 'units', char('g kg-1'));
netcdf.putAtt(nc_new, var_salinity , 'FillValue_', -1e20);
%% Temperature
netcdf.putAtt(nc_new, var_temp , 'long_name', char('sea_water_conservative_temperature'));
netcdf.putAtt(nc_new, var_temp , 'units', char('degrees celsius [c]'));
netcdf.putAtt(nc_new, var_temp , 'FillValue_', -10e20);



%% General information
str = ['Standar Atlantis tempora file - created based on the roms model and the bgm ' ...
       'configuration file. Properties averaged from hydrodynamic model outputs.'];
his = ['Created Using Javier Porobic (more details : https://github.com/jporobicg) codes on' date];
netcdf.putAtt(nc_new, NC_GLOBAL, 'Title', char(' Roms files Salish sea ocean output'));
netcdf.putAtt(nc_new, NC_GLOBAL, 'Comments', str);
netcdf.putAtt(nc_new, NC_GLOBAL, 'history', str);

% % To
%mantain control and record
%his        = ['Created by Bec Gorton and modified by Javier Porobic' date];
%nc.history = his;
salt = permute(salt, [3, 1, 2]);
temp = permute(temp, [3, 1, 2]);
w    = permute(w, [3, 1, 2]);


%% Getting sure that we are dealing with real values
w(w > 1e15)       = 0;
salt(salt > 1e15) = 0;
temp(temp > 1e15) = 0;
%% Leaving define mode
netcdf.endDef(nc_new)

%% WRITE VARIABLES
netcdf.putVar(nc_new, time, 0, size(tims, 1), tims);
netcdf.putVar(nc_new, boxes, bid);
netcdf.putVar(nc_new, level, 1 : nlv);
netcdf.putVar(nc_new, var_salinity, salt);
netcdf.putVar(nc_new, var_temp, temp);
netcdf.putVar(nc_new, var_w, w);
%% CLOSE FILES
netcdf.close(nc_new);



% %% Data
% for ii = 1:ntm
%     nc{'temperature'}(ii,:,:)  = squeeze(temp(:, ii, :));
%     nc{'salinity'}(ii,:,:)     = squeeze(salt(:, ii, :));
%     nc{'verticalflux'}(ii,:,:) = squeeze(w(:, ii, :));
% end

% nc{'boxes'}(:) = bid;
% nc{'level'}(:) = 1 : nlv;
% nc{'time'}(:)  = tims;

% %% close NETCDF file
% nc             = close(nc);

%------------------------------------------------------------------------------

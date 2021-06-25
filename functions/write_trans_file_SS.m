%  tims    Base [1990 1 1 +10], so tim = tims - greg2time([1990 1 1 10 0 0]);
%
%  fcid    If faces not just 0:(nfc-1), eg if indexed by irealfaces, then
%          fcid = (irealfaces-1)
%  fnm     output file name, incl '.nc'.  [default 'trans.nc']
%
% USAGE: write_trans_file(pt1,pt2,lr,tims,T,fcid,fnm)

function write_trans_file_SS(pt1,pt2,lr,tims,T,fcid,fnm)

nfc = size(pt1,1);

if size(T,1) ~= nfc
   error
end

%disp('Check "time" units are appropriate')

%tims
% length(tims);
% if min(tims)>10000
%     %disp('*** Converting from time0 of 1900 to 1990')
%   t90 = greg2time([1990 1 1 10 0 0]);
%   tims = tims-t90;
% end

ntm = length(tims);
if size(T,2) ~= ntm
   size(T,2)
    ntm
   error
end
nlv = size(T,3);

if nargin<6 | isempty(fcid)
   fcid = 0:(nfc-1);
end

if nargin<7 | isempty(fnm)
   fnm = 'trans.nc';
end

%% Creating file
nc_new = netcdf.create(fnm, '64BIT_OFFSET');
% --- Create global attributes
% Define dimensions:
%time      =  netcdf.defDim(nc_new, 'time', ntm);
time      =  netcdf.defDim(nc_new, 'time', netcdf.getConstant('NC_UNLIMITED'));
level     =  netcdf.defDim(nc_new, 'level', nlv);
faces     =  netcdf.defDim(nc_new, 'faces', nfc);
NC_GLOBAL = netcdf.getConstant('NC_GLOBAL'); % constant

% ncdim('time', 0, nc);
% ncdim('level', nlv, nc);
% ncdim('faces', nfc, nc);

%% DEfning variables
%% Variables and attributes
time        =  netcdf.defVar(nc_new, 'time', 'double', time);
level       =  netcdf.defVar(nc_new, 'level', 'int', level);
faces       =  netcdf.defVar(nc_new, 'faces', 'int', faces);
var_pt1x    =  netcdf.defVar(nc_new, 'pt1_x', 'float', faces);
var_pt1y    =  netcdf.defVar(nc_new, 'pt1_y', 'float', faces);
var_pt2x    =  netcdf.defVar(nc_new, 'pt2_x', 'float', faces);
var_pt2y    =  netcdf.defVar(nc_new, 'pt2_y', 'float', faces);
var_dest    =  netcdf.defVar(nc_new, 'dest_boxid', 'int', [faces]);
var_sour    =  netcdf.defVar(nc_new, 'source_boxid', 'int', [faces]);
var_transp  =  netcdf.defVar(nc_new, 'transport', 'float', [level, faces, time]);
%%var_transp  =  netcdf.defVar(nc_new, 'transport', 'float', [time, faces, level]);
%% information of variables
netcdf.putAtt(nc_new, time , 'long_name', char('simulation_time'));
netcdf.putAtt(nc_new, time , 'units', char('seconds since 1970-01-01 00:00:00'));
netcdf.putAtt(nc_new, time , 'grads_step', char('seconds'));
netcdf.putAtt(nc_new, time , 'dt', 21600);

netcdf.putAtt(nc_new, faces, 'long_name', char('face IDs'));
netcdf.putAtt(nc_new, level, 'long_name', char('layer index; 1=near-surface'));
%% pt1 x - axis
netcdf.putAtt(nc_new, var_pt1x , 'long_name', char('x coordinate of point 1 of face'));
netcdf.putAtt(nc_new, var_pt1x , 'units', char('degree_east'));
%% pt1 y - axis
netcdf.putAtt(nc_new, var_pt1y , 'long_name', char('y coordinate of point 1 of face'));
netcdf.putAtt(nc_new, var_pt1y , 'units', char('degree_north'));
%% pt2 x - axis
netcdf.putAtt(nc_new, var_pt2x , 'long_name', char('x coordinate of point 2 of face'));
netcdf.putAtt(nc_new, var_pt2x , 'units', char('degree_east'));
%% pt2 y - axis
netcdf.putAtt(nc_new, var_pt2y , 'long_name', char('y coordinate of point 2 of face'));
netcdf.putAtt(nc_new, var_pt2y , 'units', char('degree_north'));
%% Destination box
netcdf.putAtt(nc_new, var_dest , 'long_name', char('ID of destination box'));
netcdf.putAtt(nc_new, var_dest , 'units', char('IDs'));
%% Source box
netcdf.putAtt(nc_new, var_sour , 'long_name', char('ID of source box'));
netcdf.putAtt(nc_new, var_sour , 'units', char('IDs'));
%% W - velocity component
netcdf.putAtt(nc_new, var_transp , 'long_name', char('flux across face'));
netcdf.putAtt(nc_new, var_transp , 'units', char('meter second-1'));
netcdf.putAtt(nc_new, var_transp , 'comment', char('m/s is to left, viewing from pt1 to pt2'));
netcdf.putAtt(nc_new, var_transp , 'FillValue_', -1e20);
%% Documenting the file
%% General information
str = ['Transports across faces of box model derived from hydrodynamic' ...
       ' model current field. Positive numbers are water gain in dest' ...
       ' box, loss in source box.'];
his = ['Created Using Javier Porobic (more details : https://github.com/jporobicg) codes on' date];
netcdf.putAtt(nc_new, NC_GLOBAL, 'Title', char('Transport from NEMO salish sea model'));
netcdf.putAtt(nc_new, NC_GLOBAL, 'Comments', str);
netcdf.putAtt(nc_new, NC_GLOBAL, 'Institution', char('CSIRO'));
netcdf.putAtt(nc_new, NC_GLOBAL, 'Conventions', char('Standar'));
netcdf.putAtt(nc_new, NC_GLOBAL, 'history', str);

%% permutation of the dimention
T = permute(T, [3, 1, 2]);

%% getting sure that we are dealing with normal values
T(T > 1e15) = 0;
%% Leaving define mode
netcdf.endDef(nc_new)

%% WRITE VARIABLES
netcdf.putVar(nc_new, time, 0, size(tims, 1), tims);
netcdf.putVar(nc_new, level, 1 : nlv);
netcdf.putVar(nc_new, faces, fcid);
netcdf.putVar(nc_new, var_pt1x, pt1( : , 1));
netcdf.putVar(nc_new, var_pt1y, pt1( : , 2));
netcdf.putVar(nc_new, var_pt2x, pt2( : , 1));
netcdf.putVar(nc_new, var_pt2y, pt2( : , 2));
netcdf.putVar(nc_new, var_dest, lr( : , 1));
netcdf.putVar(nc_new, var_sour, lr( : , 2));
netcdf.putVar(nc_new, var_transp, T);
%% CLOSE FILES
netcdf.close(nc_new);




% nc.history = his;

% for ii = 1:ntm
%    nc{'transport'}(ii,:,:) = squeeze(T(:,ii,:));
% end

% nc{'faces'}(:) = fcid;
% nc{'level'}(:) = 1:nlv;
% nc{'pt1_x'}(:) = pt1(:,1);
% nc{'pt1_y'}(:) = pt1(:,2);
% nc{'pt2_x'}(:) = pt2(:,1);
% nc{'pt2_y'}(:) = pt2(:,2);
% nc{'dest_boxid'}(:) = lr(:,1);
% nc{'source_boxid'}(:) = lr(:,2);
% nc{'time'}(:) = tims;

% nc = close(nc);




% nc{'transport'} = ncfloat('time','faces','level');
% nc{'transport'}.long_name = 'flux across face';
% nc{'transport'}.units = 'Sverdrup';
% nc{'transport'}.comment = '+ve is to left, viewing from pt1 to pt2';
% nc{'transport'}.FillValue_ = -10e20;

% %nc('time') = 0;
% %nc('level') = nlv;
% %nc('faces') = nfc;

% nc{'time'} = ncdouble('time');
% nc{'time'}.long_name = 'time';
% nc{'time'}.units = 'days since 1990-01-01 00:00:00 +10';
% nc{'time'}.dt = 1;

% nc{'faces'} = ncint('faces');
% nc{'faces'}.long_name = 'Face IDs';

% nc{'level'} = ncint('level');
% nc{'level'}.long_name = 'layer index; 1=near-surface';
% nc{'level'}.positive = 'down';

% % other reference variables
% nc{'pt1_x'} = ncfloat('faces');
% nc{'pt1_x'}.long_name = 'x coordinate of point 1 of face';
% nc{'pt1_x'}.units = 'degree_east';

% nc{'pt1_y'} = ncfloat('faces');
% nc{'pt1_y'}.long_name = 'y coordinate of point 1 of face';
% nc{'pt1_y'}.units = 'degree_north';

% nc{'pt2_x'} = ncfloat('faces');
% nc{'pt2_x'}.long_name = 'x coordinate of point 2 of face';
% nc{'pt2_x'}.units = 'degree_east';

% nc{'pt2_y'} = ncfloat('faces');
% nc{'pt2_y'}.long_name = 'y coordinate of point 2 of face';
% nc{'pt2_y'}.units = 'degree_north';

% nc{'dest_boxid'} = ncint('faces');
% nc{'dest_boxid'}.long_name = 'ID of destination box';
% nc{'dest_boxid'}.units = 'id';

% nc{'source_boxid'} = ncint('faces');
% nc{'source_boxid'}.long_name = 'ID of source box';
% nc{'source_boxid'}.units = 'id';


% % Transports
% nc{'transport'} = ncfloat('time','faces','level');
% nc{'transport'}.long_name = 'flux across face';
% nc{'transport'}.units = 'Sverdrup';
% nc{'transport'}.comment = '+ve is to left, viewing from pt1 to pt2';
% nc{'transport'}.FillValue_ = -10e20;





% nc.Conventions = 'Standard';
% nc.institution = 'CSIRO - IMAS';
% nc.references = ' ';

% str = ['Transports across faces of box model derived from hydrodynamic' ...
%        ' model current field. Positive numbers are water gain in dest' ...
%        ' box, loss in source box.'];
% %disp(str);
% %sad = input('Type any text to add to this "title" : ','s');
% %if ~isempty(sad)
% %   str = [str sad];
% %end
% nc.title = str;
% nc.source = '';%input('Type "source" attribute: ','s');
% nc.comment = ''; %input('Type "comment" attribute: ','s');

% his = ['Created by Bec Gorton, modified by Javier Porobic, created on' date];
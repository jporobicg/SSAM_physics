%% Getting the main fluxes and variables for the SAM model  %%
%  Code created by Javier Porobic based on Jeff Dunn approach
clear
close all
addpath('/functions/')

%%%       GLOBAL VARIABLES    %%%
% Get the boxes  %%
% Read the BMG file to get the value of the parameter from the polygons
BGM_JFR_ll = '20190812_Salish_Sea_ll_fixed.bgm';
[nbox,nface,bid,cent,b_area,vert, iface, botz] = read_boxes(BGM_JFR_ll);
[nulr,nupt1,nupt2] = read_faces2(nbox, nface, bid, vert, iface, BGM_JFR_ll);
iface      = iface;  %% Id of the faces
lr         = nulr;   %% Neightbourn Layers
pt1        = nupt1;  %% Face 1
pt2        = nupt2;  %% Face 2
irealfaces = find(~isnan(nupt1(:,1)));  % (ie those ref'd in box definitions)
fcid       = (irealfaces-1);
rimn       = 10;    % default 3 is probably too few
dinc       = 1;      %% 0.1; default 10km is probably ok, but for large boxes
                      % May want to reduce the face integration step 'dinc' for
                      % models with small or narrow boxes. Because the small boxes
                      % don't have problems of hiperdifusion
dlev = [0  25 50 100 250 400 700]; %% This structure is related with the biology
                                   %% and with the maximum deph in the BMG model


%% Running the model - saving by years %%
% Transport between layers
direc = (['/datasets/work/oa-gladstone-mr/work/oceano_data_Salish_sea/Raw_transport_data/']);
%direc = (['/home/por07g/Documents/2019/Oil_spill/Hidro_Data/first_hidro/']);
files = dir([direc, 'MERGED_*.nc']); %% U and V in teh same file
cd (['/datasets/work/oa-gladstone-mr/work/oceano_data_Salish_sea/Temp']);
%% for year = 2015 : ????

fll = '/datasets/work/oa-gladstone-mr/work/oceano_data_Salish_sea/Codes_transport/mesh.nc';

for f  =  1 : length(files)
        fnm         = [direc, files(f).name];
        transport_SS(vert, pt1, pt2, dlev, dinc, rimn, fnm, fll, f);
end
look = dir(['*SS_second_Step.mat']);
for f = 1 : length(look)
    load(look(f).name)
    if f == 1
        Tfinal = T;
        nctime = tims;
    else
        Tfinal = cat(2, Tfinal, T);
        nctime = cat(1, nctime, tims);
    end
end

%% save('Tfinal.mat', 'Tfinal', 'nctime'); %% no need this bit, used just in case you
%% don??t do the transformation
guard = (['/datasets/work/oa-gladstone-mr/work/oceano_data_Salish_sea/Final/SS_Transport_2016.nc']);
write_trans_file_SS(pt1, pt2, lr, nctime, Tfinal, fcid, guard)


%% variables by layer
%% Variables
%% for y = 2015:????
varn = {'wVelocity';  'salinity';  'temperature'}
direc = ('/datasets/work/oa-gladstone-mr/work/oceano_data_Salish_sea/Raw_Variables_data/');
fll = '/datasets/work/oa-gladstone-mr/work/oceano_data_Salish_sea/Codes_transport/mesh.nc';
cd (['/datasets/work/oa-gladstone-mr/work/oceano_data_Salish_sea/Temp2/']); %% temporal foldet to store the temporal files
for v  =  1 : length(varn)
    avname  = char(varn(v));
    if~(v == 1)
        files = dir([direc, '2016-*_Raw_variables.nc']);
    else
        files = dir([direc, '2016-*_WRaw_variables.nc']);
    end

    for nfile = 1 : length(files)
        fnm   = [direc, files(nfile).name];
        box_av_SS(vert, avname, dlev, fnm, fll, nfile)
    end
    t_files = dir(['*', avname, '_SS_Second_step.mat']);
    for f = 1 : length(t_files)
        load(t_files(f).name)
        if f == 1
            Av_final = Var_avg;
            nctime   = tims;
        else
            Av_final = cat(2, Av_final, Var_avg);
            nctime   = cat(1, nctime, tims);
        end
    end
    file.save = (['Av_', avname, '.mat'])
    save(file.save, 'Av_final', 'nctime')
    % writing the NETCDF file
    %write_av_var(nctime, bid, avname, Av_final, guard)
end
%% end
% put all the variables together
temp = load(['/datasets/work/oa-gladstone-mr/work/oceano_data_Salish_sea/Temp2/Av_temperature.mat']);
salt = load(['/datasets/work/oa-gladstone-mr/work/oceano_data_Salish_sea/Temp2/Av_salinity.mat']);
vert = load(['/datasets/work/oa-gladstone-mr/work/oceano_data_Salish_sea/Temp2/Av_wVelocity.mat']);

%% Writing variables
temperature = temp.Av_final;
salinity    = salt.Av_final;
vertical    = vert.Av_final;
%% Saving the netdf file
nctime      = temp.nctime;
guard = (['/datasets/work/oa-gladstone-mr/work/oceano_data_Salish_sea/Final/SS_Variables_2016.nc']);
write_av_var_new(nctime, bid, temperature, salinity,  vertical, guard);

%% Getting the main fluxes and variables for the Salish model  %%
%  Modified by Javier Porobic
clear
close all
addpath('/datasets/work/oa-gladstone-mr/work/oceano_data_Salish_sea/Codes_transport/functions/')

%%%       GLOBAL VARIABLES    %%%
% Get the boxes  %%
% Read the BMG file to get the value of the parameter from the polygons
BGM_JFR_ll = '/datasets/work/oa-gladstone-mr/work/oceano_data_Salish_sea/Codes_transport/20190812_Salish_Sea_ll_fixed.bgm';
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
dlev = [0  25 50 100 250 400 700]; %% This structure is related with
                                             % the biology and with the
                                             % maximum deph in the BMG model

netcdf.open(fnm);
%% variables by layer
%% Variables
varn = {'microzooplankton';  'dissolved_organic_nitrogen';  'mesozooplankton'; 'diatoms'; 'ammonium'}
%%varn = {'salinity';  'temperature'}
fnm = ('/datasets/work/oa-gladstone-mr/work/oceano_data_Salish_sea/Codes_transport/variables/ubcSSg3DBiologyFields1hV19-05_3b0e_7537_e099.nc');
fll = '/datasets/work/oa-gladstone-mr/work/oceano_data_Salish_sea/Codes_transport/mesh.nc';
cd (['/datasets/work/oa-gladstone-mr/work/oceano_data_Salish_sea/Codes_transport/variables/step2']);
for v  =  1 : length(varn)
    avname  = char(varn(v));
    box_av_SS(vert, avname, dlev, fnm, fll, 1);
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






% put all the variables together

NH4 = load(['/datasets/work/oa-gladstone-mr/work/oceano_data_Salish_sea/Codes_transport/variables/step2/Av_ammonium.mat']);
DIT = load(['/datasets/work/oa-gladstone-mr/work/oceano_data_Salish_sea/Codes_transport/variables/step2/Av_diatoms.mat']);
DON = load(['/datasets/work/oa-gladstone-mr/work/oceano_data_Salish_sea/Codes_transport/variables/step2/Av_dissolved_organic_nitrogen.mat']);
MZO = load(['/datasets/work/oa-gladstone-mr/work/oceano_data_Salish_sea/Codes_transport/variables/step2/Av_mesozooplankton.mat']);
SMZ = load(['/datasets/work/oa-gladstone-mr/work/oceano_data_Salish_sea/Codes_transport/variables/step2/Av_microzooplankton.mat']);



%% Writing variables
nctime = NH4.nctime;
csvwrite('NH4.csv', squeeze(mean(NH4.Av_final, 2))); 
csvwrite('DIT.csv', squeeze(mean(DIT.Av_final, 2)));
csvwrite('DON.csv', squeeze(mean(DON.Av_final, 2)));
csvwrite('MZO.csv', squeeze(mean(MZO.Av_final, 2)));
csvwrite('SMZ.csv', squeeze(mean(SMZ.Av_final, 2)));






%% Second Set of Variables
varn = {'biogenic_silicon';  'nitrate';  'silicon'}
%%varn = {'salinity';  'temperature'}
fnm = ('/datasets/work/oa-gladstone-mr/work/oceano_data_Salish_sea/Codes_transport/variables/ubcSSg3DBiologyFields1hV19-05_47d3_3b8f_8288.nc');
fll = '/datasets/work/oa-gladstone-mr/work/oceano_data_Salish_sea/Codes_transport/mesh.nc';
cd (['/datasets/work/oa-gladstone-mr/work/oceano_data_Salish_sea/Codes_transport/variables/step2']);
for v  =  1 : length(varn)
    avname  = char(varn(v));
    box_av_SS(vert, avname, dlev, fnm, fll, 1);
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


% put all the variables together

BSI = load(['/datasets/work/oa-gladstone-mr/work/oceano_data_Salish_sea/Codes_transport/variables/step2/Av_biogenic_silicon.mat']);
NIT = load(['/datasets/work/oa-gladstone-mr/work/oceano_data_Salish_sea/Codes_transport/variables/step2/Av_nitrate.mat']);
SIL = load(['/datasets/work/oa-gladstone-mr/work/oceano_data_Salish_sea/Codes_transport/variables/step2/Av_silicon.mat']);




%% Writing variables
nctime = NH4.nctime;
csvwrite('BSI.csv', squeeze(mean(BSI.Av_final, 2))); 
csvwrite('NIT.csv', squeeze(mean(NIT.Av_final, 2)));
csvwrite('SIL.csv', squeeze(mean(SIL.Av_final, 2)));










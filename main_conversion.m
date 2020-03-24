%% Add path to needed librairies

addpath('/usr/local/share/ismrmrd/matlab/')
addpath('/home/manondesclides/Code/Valery/Code_Matlab/ismrmrd_to_nifti')
addpath('/home/manondesclides/Code/Valery/Code_Matlab/ismrmrd_to_nifti/xiangruili-dicm2nii-b76a158')
clearvars


%% IF conversion from .dat to .h5 needed (Export from : dat -> ismrmrd -> nifti)

% you should use the dedicated parameter maps for the siemens_to_ismrmrd
% conversion in order to have access to all the parameters used for nifti
% reconstruction
xmlFile = '/home/atrotier/GITHUB/ismrmrd_to_nifti/parameterMaps/IsmrmrdParameterMap_Siemens_Table.xml';
xlsFile = '/home/atrotier/GITHUB/ismrmrd_to_nifti/parameterMaps/IsmrmrdParameterMap_Siemens_Table.xsl';


[filename_list, pathname] = uigetfile('*.dat','MultiSelect','on');

idxFile = 1;
if iscell(filename_list)
    filename_tmp = filename_list{idxFile};
else
    filename_tmp = filename_list;
end

filename_raw = fullfile(pathname,filename_tmp);
[pathstr, name, ext] = fileparts(filename_raw);
% cmdStr=['siemens_to_ismrmrd -z 2 -f ' filename_raw ' -o ' pathstr '/' name '.h5'] % if my pull request is accepted
cmdStr=['siemens_to_ismrmrd -z 2 --user-map ' xmlFile ' --user-stylesheet ' xlsFile ' -f ' filename_raw ' -o ' pathstr '/' name '.h5  '];
system(cmdStr);
filename = [pathstr '/' name '.h5'];


%% ELSE : Select .h5 file (.dat is not required if you do the conversion with our parameter cards too)

%CORONAL EXAMPLES :
%R2L
% filename='/home/manondesclides/Code/Valery/Dataset5/FID/meas_MID00035_FID13784_gre3D_2_2_coronal_R2L.h5'
%F2H
%filename='/home/manondesclides/Code/Valery/Dataset5/FID/meas_MID00036_FID13785_gre3D_2_2_coronal_F2H.h5'
%F2H with rotation 
% filename='/home/manondesclides/Code/Valery/Dataset5/FID/meas_MID00045_FID13794_gre3D_2_2_coronal_F2H_rot66_72.h5'
%R2L with rotation
% filename='/home/manondesclides/Code/Valery/Dataset5/FID/meas_MID00044_FID13793_gre3D_2_2_coronal_R2L_rot32_02.h5'

% TRANVERSAL EXAMPLES : 
%A2P
%filename = '/home/manondesclides/Code/Valery/Dataset5/FID/meas_MID00033_FID13782_gre3D_2_2_tranversal_A2P.h5'
%R2L
% filename = '/home/manondesclides/Code/Valery/Dataset5/FID/meas_MID00034_FID13783_gre3D_2_2_tranversal_R2L.h5'
%A2P with rotation
%  filename = '/home/manondesclides/Code/Valery/Dataset5/FID/meas_MID00037_FID13786_gre3D_2_4_tranversal_A2P_rot11_14.h5'
%R2L with rotation
%  filename = '/home/manondesclides/Code/Valery/Dataset5/FID/meas_MID00038_FID13787_gre3D_2_2_tranversal_R2L_rot99_69.h5'
 
% SAGITTAL EXAMPLES : 
%A2P
% filename = '/home/manondesclides/Code/Valery/Dataset5/FID/meas_MID00046_FID13795_gre3D_2_4_sagittal_A2P.h5'
%H2F
%filename = '/home/manondesclides/Code/Valery/Dataset5/FID/meas_MID00047_FID13796_gre3D_2_2_sagittal_H2F.h5'
%A2P with rotation
%  filename = '/home/manondesclides/Code/Valery/Dataset5/FID/meas_MID00041_FID13790_gre3D_2_4_sagittal_A2P_rot14_05.h5'
%H2F with rotation
% filename = '/home/manondesclides/Code/Valery/Dataset5/FID/meas_MID00042_FID13791_gre3D_2_2_sagittal_H2F_rot101_38.h5'

%% Load selected file
%get filename 
[pathstr, name, ext] = fileparts(filename);
%load .h5 file and return dataset header, ismrmrd header and reconstructed
%images
[head, hdr, img_scaled] = load_ismrmrd_ifft3d_reconstruction(filename);


%% DEPRECATED : table position is stored in userParameters : table_position_offset = hdr.userParameters.userParameterLong(65).value

% Get table_position_offset : really important ! get value of the "0" position of the table in our case (SIEMENS AREA).. to complete and confirm for other systems

%need GlobalTablePosTra (table position offset) contained by siemens raw data

% [~,GlobalTablePosTra] =  system(['grep -a "GlobalTablePosTra" ' filename_raw])
% GlobalTablePosTra = regexp(GlobalTablePosTra,'<ParamLong."GlobalTablePosTra">  { [+-]?\d* + }','Match')
% GlobalTablePosTra = regexp(GlobalTablePosTra,'[+-]?\d*','Match')
% 
% table_position_offset=str2num(cell2mat(GlobalTablePosTra{1,1}));

%% Create parameters for set_nii_hdr_manon et xform_mat_manon
%needs to be done for each image (here size(img_scaled) = 1 so we work only
%on the first and only image of the dataset (in temporal dimension))

h = extract_ismrmrd_parameters_from_headers(head, hdr);



%% Get crop image, flip and rotationate to match with true Nifti image
img = flip_image(img_scaled{1});


%% Convert and save our new image in Nifti format

% create nii struct based on img
nii_empty = nii_tool('init', img); 

%in case there will be more that one temporal image in dataset :  h2 = h, with
%h a structure of size size(img_scaled)
%here, only one temporal image so h2{1} only
h2{1}=h;
pf.lefthand = 0; %to include in a function ?
[nii_filled, h3] = set_nii_hdr(nii_empty, h2, pf);

%to change ? particular case ?
fmt = 1;
rst3D = (isnumeric(fmt) && fmt>3) || (ischar(fmt) && ~isempty(regexpi(fmt, '3D')));

%save image in nifti format
nii_tool('save', nii_filled, [pathstr '/' name '_ismrmrd_to_nifti_version.nii'], rst3D);
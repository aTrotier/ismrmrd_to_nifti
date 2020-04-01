%% Add path to needed librairies

addpath('/usr/local/share/ismrmrd/matlab/')
addpath('/home/manondesclides/Code/Valery/Code_Matlab/ismrmrd_to_nifti')
addpath('/home/manondesclides/Code/Valery/Code_Matlab/ismrmrd_to_nifti/xiangruili-dicm2nii-b76a158')
clearvars


%% IF conversion from .dat to .h5 needed (Export of MP2RAGE-CS_WIP siemens from : dat -> ismrmrd -> nifti)

% filename_raw ='/media/sylvain/rawData/Project-data/SIEMENS/MP2RAGE/data_santain/CS_MP2RAGE/meas_MID00221_FID97927_sparse_mp2rage_0p8iso_acc8p5.dat';
% [pathstr, name, ext] = fileparts(filename_raw);
% cmdStr=['siemens_to_ismrmrd -z 3 -f ' filename_raw ' -o ' pathstr '/' name '.h5']
% system(cmdStr);
% filename = [pathstr '/' name '.h5'];


%% ELSE : Select .h5 file (need .dat too)

%CORONAL EXAMPLES :
%R2L
%  filename='/home/manondesclides/Data/Dataset5/FID/meas_MID00035_FID13784_gre3D_2_2_coronal_R2L.h5'
%F2H
%filename='/home/manondesclides/Data/Dataset5/FID/meas_MID00036_FID13785_gre3D_2_2_coronal_F2H.h5'
%F2H with rotation 
% filename='/home/manondesclides/Data/Dataset5/FID/meas_MID00045_FID13794_gre3D_2_2_coronal_F2H_rot66_72.h5'
%R2L with rotation
% filename='/home/manondesclides/Data/Dataset5/FID/meas_MID00044_FID13793_gre3D_2_2_coronal_R2L_rot32_02.h5'

% TRANVERSAL EXAMPLES : 
%A2P
% filename = '/home/manondesclides/Data/Dataset5/FID/meas_MID00033_FID13782_gre3D_2_2_tranversal_A2P.h5'
%R2L
% filename = '/home/manondesclides/Data/Dataset5/FID/meas_MID00034_FID13783_gre3D_2_2_tranversal_R2L.h5'
%A2P with rotation
%  filename = '/home/manondesclides/Data/Dataset5/FID/meas_MID00037_FID13786_gre3D_2_4_tranversal_A2P_rot11_14.h5'
%R2L with rotation
%  filename = '/home/manondesclides/Data/Dataset5/FID/meas_MID00038_FID13787_gre3D_2_2_tranversal_R2L_rot99_69.h5'
 
% SAGITTAL EXAMPLES : 
%A2P
% filename = '/home/manondesclides/Data/Dataset5/FID/meas_MID00046_FID13795_gre3D_2_4_sagittal_A2P.h5'
%H2F
%filename = '/home/manondesclides/Data/Dataset5/FID/meas_MID00047_FID13796_gre3D_2_2_sagittal_H2F.h5'
%A2P with rotation
%  filename = '/home/manondesclides/Data/Dataset5/FID/meas_MID00041_FID13790_gre3D_2_4_sagittal_A2P_rot14_05.h5'
%H2F with rotation
% filename = '/home/manondesclides/Data/Dataset5/FID/meas_MID00042_FID13791_gre3D_2_2_sagittal_H2F_rot101_38.h5'




%sheep_ablation0
% filename  = '/home/manondesclides/Reseau/Imagerie/DICOM_DATA/2020-02-14_Imricor3/Reco_CS_Ismrmrd/42563_42087408_42087413_914_data_from_gt_matlab.mat'
% filename_raw = '/home/manondesclides/Reseau/Imagerie/DICOM_DATA/2020-02-14_Imricor3/FID/meas_MID00914_FID16018_CS_ANGIO_1_2ISO_TESTED.dat'
%sheep_ablation1
% filename  = '/home/manondesclides/Reseau/Imagerie/DICOM_DATA/2020-02-14_Imricor3/Reco_CS_Ismrmrd/42563_45937482_45937487_41_data_from_gt_matlab.mat'
% filename_raw = '/home/manondesclides/Reseau/Imagerie/DICOM_DATA/2020-02-14_Imricor3/FID/meas_MID00041_FID16236_CS_ANGIO_1_2ISO_TESTED.dat'
%sheep_ablation2_1
% filename  = '/home/manondesclides/Reseau/Imagerie/DICOM_DATA/2020-02-14_Imricor3/Reco_CS_Ismrmrd/42563_46216010_46216015_24_data_from_gt_matlab.mat'
% filename_raw = '/home/manondesclides/Reseau/Imagerie/DICOM_DATA/2020-02-14_Imricor3/FID/meas_MID00024_FID16296_CS_ANGIO_1_2ISO_TESTED.dat'
%sheep_ablation2_2
filename  = '/home/manondesclides/Data/Dataset7/Sheep_Ablation_2/10/42563_46216010_46216029_136_data_from_gt_matlab.mat'
filename_raw = '/home/manondesclides/Data/Dataset7/Sheep_Ablation_2/10/meas_MID00136_FID16393_CS_ANGIO_1_3ISO.dat'




%% Load selected file
%get filename 
[pathstr, name, ext] = fileparts(filename);
if strcmp(ext, '.h5')
    %load .h5 file and return dataset header, ismrmrd header and reconstructed
    %images
    [head, hdr, img_scaled] = load_ismrmrd_ifft3d_reconstruction(filename);
elseif strcmp(ext, '.mat')
    [head, hdr, img_scaled] = load_reconstructed_data(filename);
end

%% Get table_position_offset : really important ! get value of the "0" position of the table in our case (SIEMENS AREA).. to complete and confirm for other systems

%need GlobalTablePosTra (table position offset) contained by siemens raw
%data


%for new version of ismrmrd (with lGlobalTablePosTra field)
if(strcmp (getfield(hdr.userParameters.userParameterDouble(2), 'name'), 'lGlobalTablePosTra'))
    table_position_offset = getfield(hdr.userParameters.userParameterDouble(2), 'value')
%for old version
else
    [~,GlobalTablePosTra] =  system(['grep -a "GlobalTablePosTra" ' filename_raw]);
    GlobalTablePosTra = regexp(GlobalTablePosTra,'<ParamLong."GlobalTablePosTra">  { [+-]?\d* + }','Match');
    GlobalTablePosTra = regexp(GlobalTablePosTra,'[+-]?\d*','Match');

    table_position_offset=str2num(cell2mat(GlobalTablePosTra{1,1}));
end


%% Create parameters for set_nii_hdr_manon et xform_mat_manon
%needs to be done for each image (here size(img_scaled) = 1 so we work only
%on the first and only image of the dataset (in temporal dimension))

h = extract_ismrmrd_parameters_from_headers(head, hdr, table_position_offset);



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
nii_tool('save', nii_filled, [pathstr '/' name '_ismrmrd_to_nifti_version_dicm2niixtrasnfor.nii'], rst3D);

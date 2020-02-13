addpath('/home/manondesclides/Code/Valery/Code_Matlab/')


%% CORONAL EXAMPLES
%R2L
% true_filename = '../Dataset5/Nifti/9_9_gre3D_2.2_coronal_R2L_20200130142618.nii';
% converted_filename = '../Dataset5/FID/meas_MID00035_FID13784_gre3D_2_2_coronal_R2L_ismrmrd_to_nifti_version.nii';
%F2H
% true_filename = '../Dataset5/Nifti/10_10_gre3D_2.2_coronal_F2H_20200130142618.nii';
% converted_filename = '../Dataset5/FID/meas_MID00036_FID13785_gre3D_2_2_coronal_F2H_ismrmrd_to_nifti_version.nii';
%F2H with rotation
% true_filename = '../Dataset5/Nifti/16_16_gre3D_2.2_coronal_F2H_rot66.72_20200130142618.nii';
% converted_filename = '../Dataset5/FID/meas_MID00045_FID13794_gre3D_2_2_coronal_F2H_rot66_72_ismrmrd_to_nifti_version.nii';
%R2L with rotation
% true_filename = '/home/manondesclides/Code/Valery/Dataset5/Nifti/15_15_gre3D_2.2_coronal_R2L_rot32.02_20200130142618.nii';
% converted_filename = '/home/manondesclides/Code/Valery/Dataset5/FID/meas_MID00044_FID13793_gre3D_2_2_coronal_R2L_rot32_02_ismrmrd_to_nifti_version.nii';


%% TRANVERSAL EXAMPLES
%A2P
% true_filename = '../Dataset5/Nifti/7_7_gre3D_2.2_tranversal_A2P_20200130142618.nii';
% converted_filename = '../Dataset5/FID/meas_MID00033_FID13782_gre3D_2_2_tranversal_A2P_ismrmrd_to_nifti_version.nii';
%R2L
% true_filename = '../Dataset5/Nifti/8_8_gre3D_2.2_tranversal_R2L_20200130142618.nii';
% converted_filename = '../Dataset5/FID/meas_MID00034_FID13783_gre3D_2_2_tranversal_R2L_ismrmrd_to_nifti_version.nii';
%A2P with rotation
% true_filename = '/home/manondesclides/Code/Valery/Dataset5/Nifti/11_11_gre3D_2.4_tranversal_A2P_rot11.14_20200130142618.nii';
% converted_filename = '/home/manondesclides/Code/Valery/Dataset5/FID/meas_MID00037_FID13786_gre3D_2_4_tranversal_A2P_rot11_14_ismrmrd_to_nifti_version.nii';
%R2L with rotation
% true_filename = '/home/manondesclides/Code/Valery/Dataset5/Nifti/12_12_gre3D_2.2_tranversal_R2L_rot99.69_20200130142618.nii';
% converted_filename = '/home/manondesclides/Code/Valery/Dataset5/FID/meas_MID00038_FID13787_gre3D_2_2_tranversal_R2L_rot99_69_ismrmrd_to_nifti_version.nii';


%% SAGITTAL EXAMPLES
%A2P
% true_filename = '../Dataset5/Nifti/17_17_gre3D_2.4_sagittal_A2P_20200130142618.nii';
% converted_filename = '../Dataset5/FID/meas_MID00046_FID13795_gre3D_2_4_sagittal_A2P_ismrmrd_to_nifti_version.nii';
%H2F
% true_filename = '../Dataset5/Nifti/18_18_gre3D_2.2_sagittal_H2F_20200130142618.nii';
% converted_filename = '../Dataset5/FID/meas_MID00047_FID13796_gre3D_2_2_sagittal_H2F_ismrmrd_to_nifti_version.nii';
%A2P rot 
% true_filename = '../Dataset5/Nifti/13_13_gre3D_2.4_sagittal_A2P_rot14.05_20200130142618.nii';
% converted_filename = '../Dataset5/FID/meas_MID00041_FID13790_gre3D_2_4_sagittal_A2P_rot14_05_ismrmrd_to_nifti_version.nii';
%H2F rot
% true_filename = '/home/manondesclides/Code/Valery/Dataset5/Nifti/14_14_gre3D_2.2_sagittal_H2F_rot101.38_20200130142618.nii';
% converted_filename = '/home/manondesclides/Code/Valery/Dataset5/FID/meas_MID00042_FID13791_gre3D_2_2_sagittal_H2F_rot101_38_ismrmrd_to_nifti_version.nii';



%% TRANSFORM MATRIX COMPARISON 
%Load true data
true_file = niftiread(true_filename);
true_info = niftiinfo(true_filename);

%Load converted data
converted_file = niftiread(converted_filename);
converted_info = niftiinfo(converted_filename);

%Print transform matrices
T_false = converted_info.Transform.T
T_true = true_info.Transform.T


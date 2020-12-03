%% Add path to needed librairies

addpath('/usr/local/share/ismrmrd/matlab/')
addpath('/.../ismrmrd_to_nifti')
addpath('/.../ismrmrd_to_nifti/xiangruili-dicm2nii-b76a158')
clearvars


%% #################### SELECT A OUPUT PATH ###############################

% If output_path = '' : default path = data_path/filename+_ismrmrd_to_nifti_version_matlab.nii

output_path = ''; % need a .nii extension



 














%% #########################################################################
% ################## DON'T CHANGE ANYTHING FROM HERE ######################
% #########################################################################


[name, pathname] = uigetfile({'*.h5;*.mat','HDF5 File, MATLAB File (ISMRMRD format) (*.h5,*.mat)'},'MultiSelect','off');
filename = fullfile(pathname,name);


% ##################### LOAD FILES AND CONVERT DATA #######################

% Load selected file
[pathstr, name, ext] = fileparts(filename);
if strcmp(ext, '.h5')
    %load .h5 file and return dataset header, ismrmrd header and reconstructed
    %images
    [head, hdr, img_scaled] = load_ismrmrd_ifft3d_reconstruction(filename);
elseif strcmp(ext, '.mat')
    %load .mat file and return dataset header, ismrmrd header and reconstructed
    %images
    [head, hdr, img_scaled] = load_reconstructed_data(filename);
end


% Create parameters for set_nii_hdr et xform_mat
% needs to be done for each image (here size(img_scaled) = 1 so we work only
% on the first and only image of the dataset (in temporal dimension))
h = extract_ismrmrd_parameters_from_headers(head, hdr);

% Get crop image, flip and rotationate to match with true Nifti image
img = ismrmrd_to_nifti.flip_image(img_scaled{1});

% Create nii struct based on img
nii_empty = nii_tool('init', img); 

% In case there will be more that one temporal image in dataset :  h2 = h, with
% h a structure of size size(img_scaled)
% here, only one temporal image so h2{1} only
h2{1}=h;
pf.lefthand = 0; %to include in a function ?
[nii_filled, h3] = ismrmrd_to_nifti.set_nii_hdr(nii_empty, h2, pf);

fmt = 1;

% Choose to save nifti dataset in one unique file or in seperate (for 4D)
rst3D = (isnumeric(fmt) && fmt>3) || (ischar(fmt) && ~isempty(regexpi(fmt, '3D')));

% Save image in nifti format
if strcmp(output_path,'')
    nii_tool('save', nii_filled, [pathstr '/' name '_ismrmrd_to_nifti_version_matlab.nii'], rst3D);
else
    nii_tool('save', nii_filled, output_path, rst3D);
end

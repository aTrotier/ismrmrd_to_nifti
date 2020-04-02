%% Add path to needed librairies

addpath('/usr/local/share/ismrmrd/matlab/')
addpath('/.../ismrmrd_to_nifti')
addpath('/.../ismrmrd_to_nifti/xiangruili-dicm2nii-b76a158')
clearvars

%% #################### SELECT A CONVERSION TYPE ###########################
% 
% ## NEW ISMRMRD VERSION
% # 0  :  IF conversion already done with the new ismrmrd version -> select the
%           .h5 or .mat
% 
% ## OLD ISMRMRD VERSION
% # 1 :   IF conversion from .dat to .h5 needed (Export from : dat -> ismrmrd -> nifti)
% #       you should use the dedicated parameter maps for the siemens_to_ismrmrd
% #       conversion in order to have access to all the parameters used for nifti
% #       reconstruction
% 
% # 2 :   IF conversion already done with our parameter maps : Select .h5
%           file or .mat
% #       (.dat is not required if you do the conversion with our parameter cards too)
% 
% # 3  :  IF conversion already done but not with our parameter maps (need .H5 or .mat and .dat)
 

conversion_type =  3;

%% #################### SELECT A OUPUT PATH ###############################

% If output_path = '' : default path = data_path/filename+_ismrmrd_to_nifti_version_matlab.nii

output_path = ''; % need a .nii extension



 














%% #########################################################################
% ################## DON'T CHANGE ANYTHING FROM HERE ######################
% #########################################################################


% ############ GIVE PATH TO DATA NEEDED FOR THE SELECTED CASE #############

% ## CASE 0, 2, 3 :
if conversion_type == 0 || conversion_type == 2 || conversion_type == 3 
    [name, pathname] = uigetfile({'*.h5;*.mat','HDF5 File, MATLAB File (ISMRMRD format) (*.h5,*.mat)'},'MultiSelect','off');
    filename = fullfile(pathname,name);
    if conversion_type == 3 
        [name_raw, pathraw] = uigetfile({'*.dat','RAW SIEMENS File (*.dat)'},'MultiSelect','off');
        filename_raw = fullfile(pathraw,name_raw)
    end
    
% ## CASE 1 :
elseif conversion_type == 1
    xmlFile = '/.../ismrmrd_to_nifti/parameterMaps/IsmrmrdParameterMap_Siemens_Table.xml'
    xlsFile = '/.../ismrmrd_to_nifti/parameterMaps/IsmrmrdParameterMap_Siemens_Table.xsl'

    [filename_list, pathname] = uigetfile({'*.dat','RAW SIEMENS File (*.dat)'},'MultiSelect','on');
    
    idxFile = 1;
    if iscell(filename_list)
        filename_tmp = filename_list{idxFile};
    else
        filename_tmp = filename_list;
    end
    
    filename_raw = fullfile(pathname,filename_tmp);
    [pathstr, name, ext] = fileparts(filename_raw);
    cmdStr=['siemens_to_ismrmrd -z 2 --user-map ' xmlFile ' --user-stylesheet ' xlsFile ' -f ' filename_raw ' -o ' pathstr '/' name '.h5  '];
    system(cmdStr);
    filename = [pathstr '/' name '.h5'];

end


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

% Correct "position" fields with the table position offset parameter
if conversion_type == 1 || conversion_type == 2
%   Find indices of GlobalTableProsTra fields
    idx_TablePos = find(strcmp({hdr.userParameters.userParameterLong.name}, 'GlobalTablePosTra')==1);
    if isempty(idx_TablePos)
        warning('Missing GlobalTablePosTra parameters in the ismrmrd dataset. Check if you use the parameter maps in the folder');
        warning("Default value for table_position_offset : 0.");
    else
        table_position_offset = hdr.userParameters.userParameterLong(idx_TablePos).value;
    end
elseif conversion_type == 3 
%   Get table_position_offset : really important ! get value of the "0" position of the table in our case (SIEMENS AREA)..
%   to complete and confirm for other systems
%   need GlobalTablePosTra (table position offset) contained by siemens raw data
    [~,GlobalTablePosTra] =  system(['grep -a "GlobalTablePosTra" ' filename_raw]);
    GlobalTablePosTra = regexp(GlobalTablePosTra,'<ParamLong."GlobalTablePosTra">  { [+-]?\d* + }','Match');
    GlobalTablePosTra = regexp(GlobalTablePosTra,'[+-]?\d*','Match');
    table_position_offset=str2num(cell2mat(GlobalTablePosTra{1,1}));
else 
    table_position_offset = 0 ;
end

% Create parameters for set_nii_hdr et xform_mat
% needs to be done for each image (here size(img_scaled) = 1 so we work only
% on the first and only image of the dataset (in temporal dimension))
h = extract_ismrmrd_parameters_from_headers(head, hdr, table_position_offset);

% Get crop image, flip and rotationate to match with true Nifti image
img = flip_image(img_scaled{1});

% Create nii struct based on img
nii_empty = nii_tool('init', img); 

% In case there will be more that one temporal image in dataset :  h2 = h, with
% h a structure of size size(img_scaled)
% here, only one temporal image so h2{1} only
h2{1}=h;
pf.lefthand = 0; %to include in a function ?
[nii_filled, h3] = set_nii_hdr(nii_empty, h2, pf);

fmt = 1;

% Choose to save nifti dataset in one unique file or in seperate (for 4D)
rst3D = (isnumeric(fmt) && fmt>3) || (ischar(fmt) && ~isempty(regexpi(fmt, '3D')));

% Save image in nifti format
if strcmp(output_path,'')
    nii_tool('save', nii_filled, [pathstr '/' name '_ismrmrd_to_nifti_version_matlab.nii'], rst3D);
else
    nii_tool('save', nii_filled, output_path, rst3D);
end

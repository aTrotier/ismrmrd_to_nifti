% data_header is the header of the rawdata. For a bbucket object -> you can find it under Bucket.data.header
% connection_header is the header in the object CONNECTION of gadgetron_matlab 
% https://github.com/gadgetron/gadgetron-matlab/blob/master/%2Bgadgetron/%2Bexternal/Connection.m

% I guess it should also work for Buffer/ReconData. Maybe take care with empty field.
function gadgetron_to_nii(img,data_header,connection_header,pathname)
% Get crop image, flip and rotationate to match with true Nifti image
img = ismrmrd_to_nifti.flip_image(abs(img));

h = extract_ismrmrd_parameters_from_headers(data_header, connection_header);

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


nii_tool('save', nii_filled, pathname, rst3D);


end
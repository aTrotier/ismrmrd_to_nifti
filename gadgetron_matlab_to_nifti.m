%
% gadgetron_matlab_to_nifti(path)
%
% Author:   Aurélien TROTIER  (a.trotier@gmail.com)
% Date:     2020-01-17
% Partner: none
% Institute: CRMSB (Bordeaux)
%
% Function description:
% Can reconstruct 2D or 3D non cartesian datasets
%
% Input:
%   img : 
%   hdr : g.xml
%   head : headers of data  extract with structure from (recon_data(1).data.headers
%
% Output:
%
%
%
%
% Algorithm & bibliography:
%
% See also :
%
% To do :
%


function gadgetron_matlab_to_nifti(img,hdr,head,s)

if isfield(s,'path')
    path=s.path;
else
    path='/tmp/';
end

save('/home/ums/Dev/ismrmrd_to_nifti/data_test/param.mat');

[ idx ] = get_non_zero_index( head )

%% create equivalent of dicom parameters in order to be used by : set_nii_hdr_aurel.m
% function are modfied version of : https://github.com/xiangruili/dicm2nii

h=struct('NumberOfTemporalPositions',1);
h.MRAcquisitionType='3D';
h.ImageOrientationPatient=[ head.phase_dir(:,idx)' -head.read_dir(:,idx)']

h.SpacingBetweenSlices=hdr.encoding.reconSpace.fieldOfView_mm.z/hdr.encoding.reconSpace.matrixSize.z;

if(hdr.encoding.reconSpace.matrixSize.z == 0)
    hdr.encoding.reconSpace.matrixSize.z = 1;
end
h.SliceThickness=hdr.encoding.reconSpace.fieldOfView_mm.z/hdr.encoding.reconSpace.matrixSize.z;

%
h.PixelSpacing = [hdr.encoding.reconSpace.fieldOfView_mm.x / hdr.encoding.reconSpace.matrixSize.x ...
    hdr.encoding.reconSpace.fieldOfView_mm.y / hdr.encoding.reconSpace.matrixSize.y ]';

% ImagePositionPatient of first image  
% Attention peut-être que faut changer les +/- 
% en fonction de HFS dans : hdr.measurementInformation.patientPosition

corner(1)=head.position(1,idx) ...
    + (hdr.encoding.reconSpace.fieldOfView_mm.x / 2.0 ) * head.read_dir(1,idx) ...
    - (hdr.encoding.reconSpace.fieldOfView_mm.y / 2.0 ) * head.phase_dir(1,idx) ...
    - (hdr.encoding.reconSpace.fieldOfView_mm.z / 2.0 ) * head.slice_dir(1,idx);

corner(2)=head.position(2,idx) ...
    + (hdr.encoding.reconSpace.fieldOfView_mm.x / 2.0 ) * head.read_dir(2,idx) ...
    - (hdr.encoding.reconSpace.fieldOfView_mm.y / 2.0 ) * head.phase_dir(2,idx) ...
    - (hdr.encoding.reconSpace.fieldOfView_mm.z / 2.0 ) * head.slice_dir(2,idx);

corner(3)=head.position(3,idx) ...
    + (hdr.encoding.reconSpace.fieldOfView_mm.x / 2.0 ) * head.read_dir(3,idx) ...
    - (hdr.encoding.reconSpace.fieldOfView_mm.y / 2.0 ) * head.phase_dir(3,idx) ...
    - (hdr.encoding.reconSpace.fieldOfView_mm.z / 2.0 ) * head.slice_dir(3,idx);

h.ImagePositionPatient = corner';

% ImagePositionPatient of last image
corner(1)=head.position(1,idx) ...
    + (hdr.encoding.reconSpace.fieldOfView_mm.x / 2.0 ) * head.read_dir(1,idx) ...
    - (hdr.encoding.reconSpace.fieldOfView_mm.y / 2.0 ) * head.phase_dir(1,idx) ...
    + (hdr.encoding.reconSpace.fieldOfView_mm.z / 2.0 ) * head.slice_dir(1,idx);

corner(2)=head.position(2,idx) ...
    + (hdr.encoding.reconSpace.fieldOfView_mm.x / 2.0 ) * head.read_dir(2,idx) ...
    - (hdr.encoding.reconSpace.fieldOfView_mm.y / 2.0 ) * head.phase_dir(2,idx) ...
    + (hdr.encoding.reconSpace.fieldOfView_mm.z / 2.0 ) * head.slice_dir(2,idx);

corner(3)=head.position(3,idx) ...
    + (hdr.encoding.reconSpace.fieldOfView_mm.x / 2.0 ) * head.read_dir(3,idx) ...
    - (hdr.encoding.reconSpace.fieldOfView_mm.y / 2.0 ) * head.phase_dir(3,idx) ...
    + (hdr.encoding.reconSpace.fieldOfView_mm.z / 2.0 ) * head.slice_dir(3,idx);

h.LastFile.ImagePositionPatient = corner';

%
h.Manufacturer=hdr.acquisitionSystemInformation.systemVendor;

%
h.Columns = hdr.encoding.reconSpace.matrixSize.y;

%
h.NiftiName=hdr.measurementInformation.protocolName;

h.InPlanePhaseEncodingDirection='COL'; %or 'ROW'

pf.lefthand = 0;

%% write nifti with header
img = permute(img, [2 1 3]);
img=flip(img,1);
img=flip(img,3);

nii = nii_tool('init', img); % create nii struct based on img

h2{1}=h;
[nii] = set_nii_hdr_aurel(nii, h2, pf);

nii.hdr.descrip=hdr.measurementInformation.measurementID;

fmt = 1;
rst3D = (isnumeric(fmt) && fmt>3) || (ischar(fmt) && ~isempty(regexpi(fmt, '3D')));
nii_tool('save', nii, [path hdr.measurementInformation.measurementID '.nii'], rst3D);


function [ output ] = get_non_zero_index( head )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    tempo=head.idx.kspace_encode_step_1;    
    idx=find(tempo~=0);      
    output=idx(1);

end

end
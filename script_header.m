%% export nifti from ismrmrd : recreate dicom entries then use a modified dicom2nifti export
load('ismrmrd_info.mat');

%% Create parameters for set_nii_hdr_aurel et xform_mat_aurel
h=struct('NumberOfTemporalPositions',1);
h.MRAcquisitionType='3D';
h.ImageOrientationPatient=[head.read_dir(:,1)' head.phase_dir(:,1)']

h.SpacingBetweenSlices=hdr.encoding.reconSpace.fieldOfView_mm.z/hdr.encoding.reconSpace.matrixSize.z;

if(hdr.encoding.reconSpace.matrixSize.z == 0)
    hdr.encoding.reconSpace.matrixSize.z = 1;
end
h.SliceThickness=hdr.encoding.reconSpace.fieldOfView_mm.z/hdr.encoding.reconSpace.matrixSize.z;

%
h.PixelSpacing = [hdr.encoding.reconSpace.fieldOfView_mm.x / hdr.encoding.reconSpace.matrixSize.x ...
    hdr.encoding.reconSpace.fieldOfView_mm.y / hdr.encoding.reconSpace.matrixSize.y ]';

%
h.Manufacturer=hdr.acquisitionSystemInformation.systemVendor;

%
h.Columns = hdr.encoding.reconSpace.matrixSize.y;

%
h.NiftiName=hdr.measurementInformation.protocolName;

h.InPlanePhaseEncodingDirection='COL'; %or 'ROW'

pf.lefthand = 0;
%% Only issue for this dataset -> I need to obtain : 
%[-103.38414726195 -165.86384794472 114.56585456728]

corner(1)=head.position(1) ...
    + (hdr.encoding.reconSpace.fieldOfView_mm.x / 2.0 ) * head.read_dir(1) ...
    - (hdr.encoding.reconSpace.fieldOfView_mm.y / 2.0 ) * head.phase_dir(1) ...
    - (hdr.encoding.reconSpace.fieldOfView_mm.z / 2.0 ) * head.slice_dir(1);

corner(2)=head.position(2) ...
    + (hdr.encoding.reconSpace.fieldOfView_mm.x / 2.0 ) * head.read_dir(2) ...
    - (hdr.encoding.reconSpace.fieldOfView_mm.y / 2.0 ) * head.phase_dir(2) ...
    - (hdr.encoding.reconSpace.fieldOfView_mm.z / 2.0 ) * head.slice_dir(2);

corner(3)=head.position(3) ...
    + (hdr.encoding.reconSpace.fieldOfView_mm.x / 2.0 ) * head.read_dir(3) ...
    - (hdr.encoding.reconSpace.fieldOfView_mm.y / 2.0 ) * head.phase_dir(3) ...
    - (hdr.encoding.reconSpace.fieldOfView_mm.z / 2.0 ) * head.slice_dir(3);

h.ImagePositionPatient = corner';


epsi = 0.5;
ImagePosition = [-103.38414726195 -165.86384794472 114.56585456728];
if (abs(h.ImagePositionPatient(1) - ImagePosition(1)) > epsi  || ...
        abs(h.ImagePositionPatient(2) - ImagePosition(2)) > epsi || ...
        abs(h.ImagePositionPatient(3) - ImagePosition(3)) > epsi)
    
    warning(['Image position calculated from ISMRMRD is different'])
    warning(['What we need : ' num2str(ImagePosition)])
    warning(['What we calculate : ' num2str(h.ImagePositionPatient')])
else
    disp('Image position patient is the same')
    warning(['What we calculate : ' num2str(h.ImagePositionPatient')])
end


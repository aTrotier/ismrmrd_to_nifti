function h = extract_ismrmrd_parameters_from_img_headers(headtmp)
%create_nifti_parameters : function that returns a structure h which
%contains parameters needed for conversion function

%% recreate head and hdr
first_index_not_empty = 1;
head.phase_dir(:,1) = headtmp.phase_dir;
head.read_dir(:,1) = headtmp.read_dir;
head.slice_dir(:,1) = headtmp.slice_dir;

hdr.encoding.reconSpace.fieldOfView_mm.x = double(headtmp.field_of_view(1));
hdr.encoding.reconSpace.fieldOfView_mm.y = double(headtmp.field_of_view(2));
hdr.encoding.reconSpace.fieldOfView_mm.z = double(headtmp.field_of_view(3));

hdr.encoding.reconSpace.matrixSize.x = double(headtmp.matrix_size(1));
hdr.encoding.reconSpace.matrixSize.y = double(headtmp.matrix_size(2));
hdr.encoding.reconSpace.matrixSize.z = double(headtmp.matrix_size(3));

head.position(:,1) = headtmp.position;


%% Create structure h
h = struct('NumberOfTemporalPositions',1);

%fill it with appropriate values
h.MRAcquisitionType = '3D';

% obscur (why a substraction ??)
h.ImageOrientationPatient = [ -head.phase_dir(:,first_index_not_empty)' -head.read_dir(:,first_index_not_empty)']

h.SpacingBetweenSlices = hdr.encoding.reconSpace.fieldOfView_mm.z/hdr.encoding.reconSpace.matrixSize.z;

if (hdr.encoding.reconSpace.matrixSize.z == 0)
    hdr.encoding.reconSpace.matrixSize.z = 1;
end

h.SliceThickness = hdr.encoding.reconSpace.fieldOfView_mm.z/hdr.encoding.reconSpace.matrixSize.z;

h.PixelSpacing = [hdr.encoding.reconSpace.fieldOfView_mm.x / hdr.encoding.reconSpace.matrixSize.x ...
                  hdr.encoding.reconSpace.fieldOfView_mm.y / hdr.encoding.reconSpace.matrixSize.y ]';
              
% % ImagePositionPatient of first slice 
% % Attention peut-être que faut changer les +/-
% % en fonction de HFS dans : hdr.measurementInformation.patientPosition


%Compute the corner position (ImagePositionPatient) of the first slice
corner_first_slice(1)= head.position(1, first_index_not_empty)  ... 
    + (hdr.encoding.reconSpace.fieldOfView_mm.x / 2.0 ) * head.read_dir(1,first_index_not_empty) ...
    + (hdr.encoding.reconSpace.fieldOfView_mm.y / 2.0 ) * head.phase_dir(1,first_index_not_empty) ...
    + (hdr.encoding.reconSpace.fieldOfView_mm.z / 2.0 ) * head.slice_dir(1,first_index_not_empty) ;

corner_first_slice(2)=head.position(2, first_index_not_empty)   ... 
    + (hdr.encoding.reconSpace.fieldOfView_mm.x / 2.0) * head.read_dir(2,first_index_not_empty) ...
    + (hdr.encoding.reconSpace.fieldOfView_mm.y / 2.0 ) * head.phase_dir(2,first_index_not_empty) ...
    + (hdr.encoding.reconSpace.fieldOfView_mm.z / 2.0) * head.slice_dir(2,first_index_not_empty) ;
   
corner_first_slice(3)=head.position(3, first_index_not_empty)  ...
    + (hdr.encoding.reconSpace.fieldOfView_mm.x / 2.0 ) * head.read_dir(3,first_index_not_empty) ...
    + (hdr.encoding.reconSpace.fieldOfView_mm.y / 2.0 ) * head.phase_dir(3,first_index_not_empty) ...
    + (hdr.encoding.reconSpace.fieldOfView_mm.z / 2.0 ) * head.slice_dir(3,first_index_not_empty) ;

h.ImagePositionPatient = corner_first_slice';

%Compute the corner position (ImagePositionPatient) of the last slice
last_index_not_empty = 1;
corner_last_slice(1)=head.position(1, first_index_not_empty)...
    + (hdr.encoding.reconSpace.fieldOfView_mm.x / 2.0 ) * head.read_dir(1,last_index_not_empty) ...
    + (hdr.encoding.reconSpace.fieldOfView_mm.y / 2.0 ) * head.phase_dir(1,last_index_not_empty) ...
    - (hdr.encoding.reconSpace.fieldOfView_mm.z / 2.0 ) * head.slice_dir(1,last_index_not_empty) ;

corner_last_slice(2)=head.position(2, first_index_not_empty) ... 
    + (hdr.encoding.reconSpace.fieldOfView_mm.x / 2.0 ) * head.read_dir(2,last_index_not_empty) ...
    + (hdr.encoding.reconSpace.fieldOfView_mm.y / 2.0 ) * head.phase_dir(2,last_index_not_empty) ...
    - (hdr.encoding.reconSpace.fieldOfView_mm.z / 2.0 ) * head.slice_dir(2,last_index_not_empty) ;

corner_last_slice(3)=head.position(3, first_index_not_empty) ...
    + (hdr.encoding.reconSpace.fieldOfView_mm.x / 2.0 ) * head.read_dir(3,last_index_not_empty) ...
    + (hdr.encoding.reconSpace.fieldOfView_mm.y / 2.0 ) * head.phase_dir(3,last_index_not_empty) ...
    - (hdr.encoding.reconSpace.fieldOfView_mm.z / 2.0 ) * head.slice_dir(3,last_index_not_empty) ;

h.LastFile.ImagePositionPatient = corner_last_slice';

h.Columns = hdr.encoding.reconSpace.matrixSize.y;

h.InPlanePhaseEncodingDirection='COL'; %or 'ROW'

%% fill weird info
h.Manufacturer='MRI Default';
h.NiftiName=num2str(headtmp.measurement_uid);
end


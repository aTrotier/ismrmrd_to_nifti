%% Export of MP2RAGE-CS_WIP siemens from : dat -> ismrmrd -> nifti 
% information about the acquistion -> set are used to store TI1 and TI2
% acquistion

%filename='/media/sylvain/rawData/Project-data/SIEMENS/MP2RAGE/data_santain/CS_MP2RAGE/meas_MID00221_FID97927_sparse_mp2rage_0p8iso_acc8p5.dat';
filename='/media/sylvain/rawData/Project-data/SIEMENS/MP2RAGE/data_santain/CS_MP2RAGE/meas_MID00222_FID97928_sparse_mp2rage_1p0iso_acc8p0.dat';
[pathstr, name, ext] = fileparts(filename);

cmdStr=['siemens_to_ismrmrd -z 3 -f ' filename ' -o ' pathstr '/' name '.h5']
system(cmdStr);


filename = [pathstr '/' name '.h5'];


if exist(filename, 'file')
    dset = ismrmrd.Dataset(filename, 'dataset');
else
    error(['File ' filename ' does not exist.  Please generate it.'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read some fields from the XML header %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We need to check if optional fields exists before trying to read them

hdr = ismrmrd.xml.deserialize(dset.readxml);

D = dset.readAcquisition();
head=D.head;
%% Encoding and reconstruction information
% Matrix size
enc_Nx = hdr.encoding.encodedSpace.matrixSize.x;
enc_Ny = hdr.encoding.encodedSpace.matrixSize.y;
enc_Nz = hdr.encoding.encodedSpace.matrixSize.z;
rec_Nx = hdr.encoding.reconSpace.matrixSize.x;
rec_Ny = hdr.encoding.reconSpace.matrixSize.y;
rec_Nz = hdr.encoding.reconSpace.matrixSize.z;

% Field of View
enc_FOVx = hdr.encoding.encodedSpace.fieldOfView_mm.x;
enc_FOVy = hdr.encoding.encodedSpace.fieldOfView_mm.y;
enc_FOVz = hdr.encoding.encodedSpace.fieldOfView_mm.z;
rec_FOVx = hdr.encoding.reconSpace.fieldOfView_mm.x;
rec_FOVy = hdr.encoding.reconSpace.fieldOfView_mm.y;
rec_FOVz = hdr.encoding.reconSpace.fieldOfView_mm.z;

% Number of slices, coils, repetitions, contrasts etc.
% We have to wrap the following in a try/catch because a valid xml header may
% not have an entry for some of the parameters

try
    nSlices = hdr.encoding.encodingLimits.slice.maximum + 1;
catch
    nSlices = 1;
end

try
    nCoils = hdr.acquisitionSystemInformation.receiverChannels;
catch
    nCoils = 1;
end

try
    nReps = hdr.encoding.encodingLimits.repetition.maximum + 1;
catch
    nReps = 1;
end

try
    nContrasts = hdr.encoding.encodingLimits.contrast.maximum + 1 + 1;
catch
    nContrasts = 1;
end

try
    set = hdr.encoding.encodingLimits.set +1;
catch
    set = 1;
end
%% Ignore noise scans
% TODO add a pre-whitening example
% Find the first non-noise scan
% This is how to check if a flag is set in the acquisition header
isNoise = D.head.flagIsSet('ACQ_IS_NOISE_MEASUREMENT');
firstScan = find(isNoise==0,1,'first');
if firstScan > 1
    noise = D.select(1:firstScan-1);
else
    noise = [];
end
meas  = D.select(firstScan:D.getNumber);

%clear D
%% Reconstruct images
% Since the entire file is in memory we can use random access
% Loop over repetitions, contrasts, slices
reconImages = {};
nimages = 0;
for rep = 1:nReps
    for nset = 1:set %store TI1 and TI2 acquisition
        for slice = 1:nSlices
            % Initialize the K-space storage array
            K = zeros(enc_Nx, enc_Ny, enc_Nz, nCoils);
            % Select the appropriate measurements from the data
            acqs = find(  (meas.head.idx.set==(nset-1)) ...
                        & (meas.head.idx.repetition==(rep-1)) ...
                        & (meas.head.idx.slice==(slice-1)) ...
                        & (~ meas.head.flagIsSet('ACQ_IS_NAVIGATION_DATA')) ...
                        & (~ meas.head.flagIsSet('ACQ_IS_PHASECORR_DATA')) ...
                        & (~ meas.head.flagIsSet('ACQ_IS_PARALLEL_CALIBRATION')) ...
                        & (~ meas.head.flagIsSet('ACQ_IS_NOISE_MEASUREMENT')) ...
                        & (~ meas.head.flagIsSet('ACQ_IS_DUMMYSCAN_DATA')) ...
                        & (~ meas.head.flagIsSet('ACQ_IS_RTFEEDBACK_DATA')) ...
                        & (~ meas.head.flagIsSet('ACQ_IS_HPFEEDBACK_DATA')) ...
                        );
            for p = 1:length(acqs)
                ky = meas.head.idx.kspace_encode_step_1(acqs(p)) + 1;
                kz = meas.head.idx.kspace_encode_step_2(acqs(p)) + 1;
                K(:,ky,kz,:) = meas.data{acqs(p)};
            end
            % Reconstruct in x
            K = fftshift(ifft(fftshift(K,1),[],1),1);
            % Chop if needed
            if (enc_Nx == rec_Nx)
                im = K;
            else
                ind1 = floor((enc_Nx - rec_Nx)/2)+1;
                ind2 = floor((enc_Nx - rec_Nx)/2)+rec_Nx;
                im = K(ind1:ind2,:,:,:);
            end
            % Reconstruct in y then z
            im = fftshift(ifft(fftshift(im,2),[],2),2);
            if size(im,3)>1
                im = fftshift(ifft(fftshift(im,3),[],3),3);
            end
            
            % Combine SOS across coils
            im = sqrt(sum(abs(im).^2,4));
            
            % Append
            nimages = nimages + 1;
            reconImages{nimages} = im;
        end
    end
end

% Display the first image
figure
colormap gray
subplot(2,1,1);
imagesc(reconImages{1}(:,:,end/2)); axis image; axis off; colorbar;
subplot(2,1,2);
imagesc(reconImages{2}(:,:,end/2)); axis image; axis off; colorbar;
%%
head.read_dir(:,1)
head.phase_dir(:,1)
head.slice_dir(:,1)


%% Create parameters for set_nii_hdr_aurel et xform_mat_aurel


h=struct('NumberOfTemporalPositions',1);
h.MRAcquisitionType='3D';
h.ImageOrientationPatient=[ head.phase_dir(:,1)' -head.read_dir(:,1)']

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

% ImagePositionPatient of last image
corner(1)=head.position(1) ...
    + (hdr.encoding.reconSpace.fieldOfView_mm.x / 2.0 ) * head.read_dir(1) ...
    - (hdr.encoding.reconSpace.fieldOfView_mm.y / 2.0 ) * head.phase_dir(1) ...
    + (hdr.encoding.reconSpace.fieldOfView_mm.z / 2.0 ) * head.slice_dir(1);

corner(2)=head.position(2) ...
    + (hdr.encoding.reconSpace.fieldOfView_mm.x / 2.0 ) * head.read_dir(2) ...
    - (hdr.encoding.reconSpace.fieldOfView_mm.y / 2.0 ) * head.phase_dir(2) ...
    + (hdr.encoding.reconSpace.fieldOfView_mm.z / 2.0 ) * head.slice_dir(2);

corner(3)=head.position(3) ...
    + (hdr.encoding.reconSpace.fieldOfView_mm.x / 2.0 ) * head.read_dir(3) ...
    - (hdr.encoding.reconSpace.fieldOfView_mm.y / 2.0 ) * head.phase_dir(3) ...
    + (hdr.encoding.reconSpace.fieldOfView_mm.z / 2.0 ) * head.slice_dir(3);

h.LastFile.ImagePositionPatient = corner';

%
h.Manufacturer=hdr.acquisitionSystemInformation.systemVendor;

%
h.Columns = hdr.encoding.reconSpace.matrixSize.y;

%
h.NiftiName=hdr.measurementInformation.protocolName;

h.InPlanePhaseEncodingDirection='COL'; %or 'ROW'

pf.lefthand = 0;

%%
img=reconImages{1};
img = permute(img, [2 1 3]);
nii = nii_tool('init', img); % create nii struct based on img

h2{1}=h;
[nii] = set_nii_hdr_aurel(nii, h2, pf);

%%
fmt = 1;

rst3D = (isnumeric(fmt) && fmt>3) || (ischar(fmt) && ~isempty(regexpi(fmt, '3D')));
nii_tool('save', nii, [pathstr '/' name '.nii'], rst3D);
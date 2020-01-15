% export in ismrmrd
filename='/media/sylvain/rawData/Project-data/SIEMENS/MP2RAGE/data_santain/CS_MP2RAGE/meas_MID00221_FID97927_sparse_mp2rage_0p8iso_acc8p5.dat';
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

clear D dset
%% Reconstruct images
% Since the entire file is in memory we can use random access
% Loop over repetitions, contrasts, slices

nContrasts=1;

reconImages = {};
nimages = 0;
for rep = 1:nReps
    for contrast = 1:nContrasts
        for slice = 1:nSlices
            % Initialize the K-space storage array
            K = single(zeros(enc_Nx, enc_Ny, enc_Nz, 1));
            % Select the appropriate measurements from the data
            acqs = find(  (meas.head.idx.contrast==(contrast-1)) ...
                & (meas.head.idx.repetition==(rep-1)) ...
                & (meas.head.idx.set==(set)) ...
                & (meas.head.idx.slice==(slice-1)));
            for p = 1:length(acqs)
                ky = meas.head.idx.kspace_encode_step_1(acqs(p)) + 1;
                kz = meas.head.idx.kspace_encode_step_2(acqs(p)) + 1;
                %K(:,ky,kz,:) = single(meas.data{acqs(p)});
                K(:,ky,kz,:) = single(meas.data{acqs(p)}(:,:,:,1));
            end
            % Reconstruct in x
            K = fftshift(ifft(fftshift(K,1),[],1),1);
            % Chop if needed
            if (enc_Nx == rec_Nx)
                kspace = K;
            else
                ind1 = floor((enc_Nx - rec_Nx)/2)+1;
                ind2 = floor((enc_Nx - rec_Nx)/2)+rec_Nx;
                kspace = K(ind1:ind2,:,:,:);
            end
            clear K meas
            % go back to kspace in x
            kspace = fftshift(fft(fftshift(kspace,1),[],1),1);
            
            
        end
    end
end

img  = fftshift(fft(fftshift(kspace,[1 2 3]),[],[1 2 3]),[1 2 3]);
%%
head.read_dir(:,1)
head.phase_dir(:,1)
head.slice_dir(:,1)


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
corner(1)=head.position(1) - ...
    (hdr.encoding.reconSpace.fieldOfView_mm.x / 2.0 ) * head.read_dir(1) - ...
    (hdr.encoding.reconSpace.fieldOfView_mm.y / 2.0 ) * head.phase_dir(1) - ...
    (hdr.encoding.reconSpace.fieldOfView_mm.z / 2.0 ) * head.slice_dir(1);

corner(2)=head.position(2) - ...
    (hdr.encoding.reconSpace.fieldOfView_mm.x / 2.0 ) * head.read_dir(2) - ...
    (hdr.encoding.reconSpace.fieldOfView_mm.y / 2.0 ) * head.phase_dir(2) - ...
    (hdr.encoding.reconSpace.fieldOfView_mm.z / 2.0 ) * head.slice_dir(2);

corner(3)=head.position(3) - ...
    (hdr.encoding.reconSpace.fieldOfView_mm.x / 2.0 ) * head.read_dir(3) - ...
    (hdr.encoding.reconSpace.fieldOfView_mm.y / 2.0 ) * head.phase_dir(3) - ...
    (hdr.encoding.reconSpace.fieldOfView_mm.z / 2.0 ) * head.slice_dir(3);

h.ImagePositionPatient = corner';
h.ImagePositionPatient

%
h.Manufacturer=hdr.acquisitionSystemInformation.systemVendor;

%
h.Columns = hdr.encoding.reconSpace.matrixSize.y;

%
h.NiftiName=hdr.measurementInformation.protocolName;

h.InPlanePhaseEncodingDirection='COL'; %or 'ROW'

pf.lefthand = 0;

%%
%img = zeros(hdr.encoding.reconSpace.matrixSize.x,hdr.encoding.reconSpace.matrixSize.y,hdr.encoding.reconSpace.matrixSize.z);
img=flip(flip(flip(T1mapCS,1),2),3);
nii = nii_tool('init', img); % create nii struct based on img

h2{1}=h;
[nii, h] = set_nii_hdr_aurel(nii, h2, pf);

%%
fmt = 1;

rst3D = (isnumeric(fmt) && fmt>3) || (ischar(fmt) && ~isempty(regexpi(fmt, '3D')));
nii_tool('save', nii, [pathstr '/' name '.nii'], rst3D);
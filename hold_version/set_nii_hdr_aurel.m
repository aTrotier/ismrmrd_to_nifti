function [nii, h] = set_nii_hdr_aurel(nii, h, pf)
dim = nii.hdr.dim(2:4); nVol = nii.hdr.dim(5);
fld = 'NumberOfTemporalPositions';
if ~isfield(h{1}, fld) && nVol>1, h{1}.(fld) = nVol; end

% Transformation matrix: most important feature for nii
[ixyz, R, pixdim, xyz_unit] = xform_mat_aurel(h{1}, dim); % R: dicom xform matrix
R(1:2,:) = -R(1:2,:); % dicom LPS to nifti RAS, xform matrix before reorient

% Store CardiacTriggerDelayTime
fld = 'CardiacTriggerDelayTime';
if ~isfield(h{1}, 'CardiacTriggerDelayTimes') && nVol>1 && isfield(h{1}, fld)
    if numel(h) == 1 % multi frames
        iFrames = 1:dim(3):dim(3)*nVol;
        if isfield(h{1}, 'SortFrames'), iFrames = h{1}.SortFrames(iFrames); end
        s2 = struct(fld, nan(1,nVol));
        s2 = dicm_hdr(h{1}, s2, iFrames);
        tt = s2.(fld);
    else
        tt = zeros(1, nVol);
        inc = numel(h) / nVol;
        for j = 1:nVol
            tt(j) = tryGetField(h{(j-1)*inc+1}, fld, 0);
        end
    end
    if ~all(diff(tt)==0), h{1}.CardiacTriggerDelayTimes = tt; end
end

% Get EchoTime for each vol
if ~isfield(h{1}, 'EchoTimes') && nVol>1 && isfield(h{1}, 'EchoTime')
    if numel(h) == 1 % 4D multi frames
        iFrames = 1:dim(3):dim(3)*nVol;
        if isfield(h{1}, 'SortFrames'), iFrames = h{1}.SortFrames(iFrames); end
        s2 = struct('EffectiveEchoTime', nan(1,nVol));
        s2 = dicm_hdr(h{1}, s2, iFrames);
        ETs = s2.EffectiveEchoTime;
    else % regular dicom. Vida done previously
        ETs = zeros(1, nVol);
        inc = numel(h) / nVol;
        for j = 1:nVol
            ETs(j) = tryGetField(h{(j-1)*inc+1}, 'EchoTime', 0);
        end
    end
    if ~all(diff(ETs)==0), h{1}.EchoTimes = ETs; end
end

% set TR and slice timing related info before re-orient
%[h, nii.hdr] = sliceTiming(h, nii.hdr);
nii.hdr.xyzt_units = xyz_unit + nii.hdr.xyzt_units; % normally: mm (2) + sec (8)
s = h{1};

% Store motion parameters for MoCo series
if all(isfield(s, {'RBMoCoTrans' 'RBMoCoRot'})) && nVol>1
    inc = numel(h) / nVol;
    s.RBMoCoTrans = zeros(nVol, 3);
    s.RBMoCoRot   = zeros(nVol, 3);
    for j = 1:nVol
        s.RBMoCoTrans(j,:) = tryGetField(h{(j-1)*inc+1}, 'RBMoCoTrans', [0 0 0]);
        s.RBMoCoRot(j,:)   = tryGetField(h{(j-1)*inc+1}, 'RBMoCoRot',   [0 0 0]);
    end
end

% Store FrameReferenceTime: seen in Philips PET
if isfield(s, 'FrameReferenceTime') && nVol>1
    inc = numel(h) / nVol;
    vTime = zeros(1, nVol);
    dict = dicm_dict('', 'FrameReferenceTime');
    for j = 1:nVol
        s2 = dicm_hdr(h{(j-1)*inc+1}.Filename, dict);
        vTime(j) = tryGetField(s2, 'FrameReferenceTime', 0);
    end
    if vTime(1) > vTime(end) % could also re-read sorted h{i}{1}
        vTime = flip(vTime);
        nii.img = flip(nii.img, 4);
    end
    s.VolumeTiming = vTime / 1000; % ms to seconds
end

% dim_info byte: freq_dim, phase_dim, slice_dim low to high, each 2 bits
[phPos, iPhase] = phaseDirection(s); % phPos relative to image in FSL feat!
if     iPhase == 2, fps_bits = [1 4 16];
elseif iPhase == 1, fps_bits = [4 1 16]; 
else,               fps_bits = [0 0 16];
end

%% Reorient if MRAcquisitionType==3D || isDTI && nSL>1
% If FSL etc can read dim_info for STC, we can always reorient.
[~, perm] = sort(ixyz); % may permute 3 dimensions in this order
if (strcmp(tryGetField(s, 'MRAcquisitionType', ''), '3D') || s.isDTI) && ...
        dim(3)>1 && (~isequal(perm, 1:3)) % skip if already XYZ order
    R(:, 1:3) = R(:, perm); % xform matrix after perm
    fps_bits = fps_bits(perm);
    ixyz = ixyz(perm); % 1:3 after perm
    dim = dim(perm);
    pixdim = pixdim(perm);
    nii.hdr.dim(2:4) = dim;
    nii.img = permute(nii.img, [perm 4:8]);
    if isfield(s, 'bvec'), s.bvec = s.bvec(:, perm); end
end
iSL = find(fps_bits==16);
iPhase = find(fps_bits==4); % axis index for phase_dim in re-oriented img

nii.hdr.dim_info = (1:3) * fps_bits'; % useful for EPI only
nii.hdr.pixdim(2:4) = pixdim; % voxel zize

flp = R(ixyz+[0 3 6])<0; % flip an axis if true
d = det(R(:,1:3)) * prod(1-flp*2); % det after all 3 axis positive
if (d>0 && pf.lefthand) || (d<0 && ~pf.lefthand)
    flp(1) = ~flp(1); % left or right storage
end
rotM = diag([1-flp*2 1]); % 1 or -1 on diagnal
rotM(1:3, 4) = (dim-1) .* flp; % 0 or dim-1
R = R / rotM; % xform matrix after flip
for k = 1:3, if flp(k), nii.img = flip(nii.img, k); end; end
if flp(iPhase), phPos = ~phPos; end
if isfield(s, 'bvec'), s.bvec(:, flp) = -s.bvec(:, flp); end
if flp(iSL) && isfield(s, 'SliceTiming') % slices flipped
    s.SliceTiming = flip(s.SliceTiming);
    sc = nii.hdr.slice_code;
    if sc>0, nii.hdr.slice_code = sc+mod(sc,2)*2-1; end % 1<->2, 3<->4, 5<->6
end

% sform
frmCode = all(isfield(s, {'ImageOrientationPatient' 'ImagePositionPatient'}));
frmCode = tryGetField(s, 'TemplateSpace', frmCode);
nii.hdr.sform_code = frmCode; % 1: SCANNER_ANAT
nii.hdr.srow_x = R(1,:);
nii.hdr.srow_y = R(2,:);
nii.hdr.srow_z = R(3,:);

R0 = normc(R(:, 1:3));
sNorm = null(R0(:, setdiff(1:3, iSL))');
if sign(sNorm(ixyz(iSL))) ~= sign(R(ixyz(iSL),iSL)), sNorm = -sNorm; end
shear = norm(R0(:,iSL)-sNorm) > 0.01;
R0(:,iSL) = sNorm;

% qform
nii.hdr.qform_code = frmCode;
nii.hdr.qoffset_x = R(1,4);
nii.hdr.qoffset_y = R(2,4);
nii.hdr.qoffset_z = R(3,4);
[q, nii.hdr.pixdim(1)] = dcm2quat(R0); % 3x3 dir cos matrix to quaternion
nii.hdr.quatern_b = q(2);
nii.hdr.quatern_c = q(3);
nii.hdr.quatern_d = q(4);

if shear
    nii.hdr.hdrTilt = nii.hdr; % copy all hdr for tilt version
    nii.hdr.qform_code = 0; % disable qform
    gantry = tryGetField(s, 'GantryDetectorTilt', 0);
    nii.hdr.hdrTilt.pixdim(iSL+1) = norm(R(1:3, iSL)) * cosd(gantry);
    R(1:3, iSL) = sNorm * nii.hdr.hdrTilt.pixdim(iSL+1);
    nii.hdr.hdrTilt.srow_x = R(1,:);
    nii.hdr.hdrTilt.srow_y = R(2,:);
    nii.hdr.hdrTilt.srow_z = R(3,:);
end

% store some possibly useful info in descrip and other text fields
str = tryGetField(s, 'ImageComments', '');
if isType(s, '\MOCO\'), str = ''; end % useless for MoCo
foo = tryGetField(s, 'StudyComments');
if ~isempty(foo), str = [str ';' foo]; end
str = [str ';' sscanf(s.Manufacturer, '%s', 1)];
foo = tryGetField(s, 'ProtocolName');
if ~isempty(foo), str = [str ';' foo]; end
nii.hdr.aux_file = str; % char[24], info only


if ~isempty(iPhase)
    if isempty(phPos), pm = '?'; b67 = 0;
    elseif phPos,      pm = '';  b67 = 1;
    else,              pm = '-'; b67 = 2;
    end
    nii.hdr.dim_info = nii.hdr.dim_info + b67*64;
    axes = 'xyz'; % actually ijk
    phDir = [pm axes(iPhase)];
    s.UnwarpDirection = phDir;

end


% slope and intercept: apply to img if no rounding error 
sclApplied = tryGetField(s, 'ApplyRescale', false);
if any(isfield(s, {'RescaleSlope' 'RescaleIntercept'})) && ~sclApplied
    slope = tryGetField(s, 'RescaleSlope', 1); 
    inter = tryGetField(s, 'RescaleIntercept', 0);
    if isfield(s, 'MRScaleSlope') % Philips: see PAR file for detail
        inter = inter / (slope * double(s.MRScaleSlope));
        slope = 1 / double(s.MRScaleSlope);
    end
    val = sort(double([max(nii.img(:)) min(nii.img(:))]) * slope + inter);
    dClass = class(nii.img);
    if isa(nii.img, 'float') || (mod(slope,1)==0 && mod(inter,1)==0 ... 
            && val(1)>=intmin(dClass) && val(2)<=intmax(dClass))
        nii.img = nii.img * slope + inter; % apply to img if no rounding
    else
        nii.hdr.scl_slope = slope;
        nii.hdr.scl_inter = inter;
    end
elseif sclApplied && isfield(s, 'MRScaleSlope')
    slope = tryGetField(s, 'RescaleSlope', 1) * s.MRScaleSlope; 
    nii.img = nii.img / slope;
end

% if pf.scale_16bit && any(nii.hdr.datatype==[4 512]) % like dcm2niix
%     if nii.hdr.datatype == 4 % int16
%         scale = floor(32000 / double(max(abs(nii.img(:)))));
%     else % datatype==512 % uint16
%         scale = floor(64000 / double((max(nii.img(:)))));
%     end
%     nii.img = nii.img * scale;
%     nii.hdr.scl_slope = nii.hdr.scl_slope / scale;
% end
h{1} = s;

function [phPos, iPhase] = phaseDirection(s)
phPos = []; iPhase = [];
fld = 'InPlanePhaseEncodingDirection';
if isfield(s, fld)
    if     strncmpi(s.(fld), 'COL', 3), iPhase = 2; % based on dicm_img(s,0)
    elseif strncmpi(s.(fld), 'ROW', 3), iPhase = 1;
    else, errorLog(['Unknown ' fld ' for ' s.NiftiName ': ' s.(fld)]);
    end
end

if isfield(s, 'CSAImageHeaderInfo') % SIEMENS
    phPos = csa_header(s, 'PhaseEncodingDirectionPositive'); % image ref
% elseif isfield(s, 'ProtocolDataBlock') % GE
%     try % VIEWORDER "1" == bottom_up
%         phPos = s.ProtocolDataBlock.VIEWORDER == '1';
%     end
elseif isfield(s, 'UserDefineData') % GE
    % https://github.com/rordenlab/dcm2niix/issues/163
    try
    b = s.UserDefineData;
    i = typecast(b(25:26), 'uint16'); % hdr_offset
    v = typecast(b(i+1:i+4), 'single'); % 5.0 to 40.0
    if v >= 25.002, i = i + 76; flag2_off = 777; else, flag2_off = 917; end
    sliceOrderFlag = bitget(b(i+flag2_off), 2);
    phasePolarFlag = bitget(b(i+49), 3);
    phPos = ~xor(phasePolarFlag, sliceOrderFlag);
    end
else
    if isfield(s, 'Stack') % Philips
        try d = s.Stack.Item_1.MRStackPreparationDirection(1); catch, return; end
    elseif isfield(s, 'PEDirectionDisplayed') % UIH
        try d = s.PEDirectionDisplayed(1); catch, return; end
    elseif isfield(s, 'Private_0177_1100') % Bruker
        expr ='(?<=\<\+?)[LRAPSI]{1}(?=;\s*phase\>)'; % <+P;phase> or <P;phase>  
        d = regexp(char(s.Private_0177_1100'), expr, 'match', 'once');
        id = regexp('LRAPSI', d);
        id = id + mod(id,2)*2-1;
        str = 'LRAPFH'; d = str(id);
    else % unknown Manufacturer
        return;
    end
    try R = reshape(s.ImageOrientationPatient, 3, 2); catch, return; end
    [~, ixy] = max(abs(R)); % like [1 2]
    if isempty(iPhase) % if no InPlanePhaseEncodingDirection
        iPhase = strfind('RLAPFH', d);
        iPhase = ceil(iPhase/2); % 1/2/3 for RL/AP/FH
        iPhase = find(ixy==iPhase); % now 1 or 2
    end
    if     any(d == 'LPH'), phPos = false; % in dicom ref
    elseif any(d == 'RAF'), phPos = true;
    end
    if R(ixy(iPhase), iPhase)<0, phPos = ~phPos; end % tricky! in image ref
end

%% Subfunction: get field if exist, return default value otherwise
function val = tryGetField(s, field, dftVal)
if isfield(s, field), val = s.(field); 
elseif nargin>2, val = dftVal;
else, val = [];
end

%% 
function v = normc(M)
vn = sqrt(sum(M .^ 2)); % vn = vecnorm(M);
vn(vn==0) = 1;
v = bsxfun(@rdivide, M, vn);

%% Subfunction, Convert 3x3 direction cosine matrix to quaternion
% Simplied from Quaternions by Przemyslaw Baranski 
function [q, proper] = dcm2quat(R)
% [q, proper] = dcm2quat(R)
% Retrun quaternion abcd from normalized matrix R (3x3)
proper = sign(det(R));
if proper<0, R(:,3) = -R(:,3); end

q = sqrt([1 1 1; 1 -1 -1; -1 1 -1; -1 -1 1] * diag(R) + 1) / 2;
if ~isreal(q(1)), q(1) = 0; end % if trace(R)+1<0, zero it
[mx, ind] = max(q);
mx = mx * 4;

if ind == 1
    q(2) = (R(3,2) - R(2,3)) /mx;
    q(3) = (R(1,3) - R(3,1)) /mx;
    q(4) = (R(2,1) - R(1,2)) /mx;
elseif ind ==  2
    q(1) = (R(3,2) - R(2,3)) /mx;
    q(3) = (R(1,2) + R(2,1)) /mx;
    q(4) = (R(3,1) + R(1,3)) /mx;
elseif ind == 3
    q(1) = (R(1,3) - R(3,1)) /mx;
    q(2) = (R(1,2) + R(2,1)) /mx;
    q(4) = (R(2,3) + R(3,2)) /mx;
elseif ind == 4
    q(1) = (R(2,1) - R(1,2)) /mx;
    q(2) = (R(3,1) + R(1,3)) /mx;
    q(3) = (R(2,3) + R(3,2)) /mx;
end
if q(1)<0, q = -q; end % as MRICron


%% Subfunction: get dicom xform matrix and related info
function [ixyz, R, pixdim, xyz_unit] = xform_mat_aurel(s, dim)
haveIOP = isfield(s, 'ImageOrientationPatient');
if haveIOP, R = reshape(s.ImageOrientationPatient, 3, 2);
else, R = [1 0 0; 0 1 0]';
end
R(:,3) = cross(R(:,1), R(:,2)); % right handed, but sign may be wrong
foo = abs(R);
[~, ixyz] = max(foo); % orientation info: perm of 1:3
if ixyz(2) == ixyz(1), foo(ixyz(2),2) = 0; [~, ixyz(2)] = max(foo(:,2)); end
if any(ixyz(3) == ixyz(1:2)), ixyz(3) = setdiff(1:3, ixyz(1:2)); end
if nargout<2, return; end
iSL = ixyz(3); % 1/2/3 for Sag/Cor/Tra slice
signSL = sign(R(iSL, 3));

try 
    pixdim = s.PixelSpacing([2 1]);
    xyz_unit = 2; % mm
catch
    pixdim = [1 1]'; % fake
    xyz_unit = 0; % no unit information
end
thk = tryGetField(s, 'SpacingBetweenSlices');
if isempty(thk), thk = tryGetField(s, 'SliceThickness', pixdim(1)); end
pixdim = [pixdim; thk];
haveIPP = isfield(s, 'ImagePositionPatient');
if haveIPP, ipp = s.ImagePositionPatient; else, ipp = -(dim'.* pixdim)/2; end
% Next is almost dicom xform matrix, except mosaic trans and unsure slice_dir
R = [R * diag(pixdim) ipp];

if dim(3)<2, return; end % don't care direction for single slice

if s.Columns>dim(1) && ~strncmpi(s.Manufacturer, 'UIH', 3) % Siemens mosaic
    R(:,4) = R * [ceil(sqrt(dim(3))-1)*dim(1:2)/2 0 1]'; % real slice location
    vec = csa_header(s, 'SliceNormalVector'); % mosaic has this
    if ~isempty(vec) % exist for all tested data
        if sign(vec(iSL)) ~= signSL, R(:,3) = -R(:,3); end
        return;
    end
elseif isfield(s, 'LastFile') && isfield(s.LastFile, 'ImagePositionPatient')
    R(:, 3) = (s.LastFile.ImagePositionPatient - R(:,4)) / (dim(3)-1);
    thk = norm(R(:,3)); % override slice thickness if it is off
    if abs(pixdim(3)-thk)/thk > 0.01, pixdim(3) = thk; warning('slice thikness does not correspond to R'); end
    return; % almost all non-mosaic images return from here
end

% Rest of the code is almost unreachable
if isfield(s, 'CSASeriesHeaderInfo') % Siemens both mosaic and regular
    ori = {'Sag' 'Cor' 'Tra'}; ori = ori{iSL};
    sNormal = asc_header(s, ['sSliceArray.asSlice[0].sNormal.d' ori]);
    if asc_header(s, ['sSliceArray.ucImageNumb' ori]), sNormal = -sNormal; end
    if sign(sNormal) ~= signSL, R(:,3) = -R(:,3); end
    if ~isempty(sNormal), return; end
end

pos = []; % volume center we try to retrieve
if isfield(s, 'LastScanLoc') && isfield(s, 'FirstScanLocation') % GE
    pos = (s.LastScanLoc + s.FirstScanLocation) / 2; % mid-slice center
    if iSL<3, pos = -pos; end % RAS convention!
    pos = pos - R(iSL, 1:2) * (dim(1:2)'-1)/2; % mid-slice location
end

if isempty(pos) && isfield(s, 'Stack') % Philips
    ori = {'RL' 'AP' 'FH'}; ori = ori{iSL};
    pos = tryGetField(s.Stack.Item_1, ['MRStackOffcentre' ori]);
    pos = pos - R(iSL, 1:2) * dim(1:2)'/2; % mid-slice location
end

if isempty(pos) % keep right-handed, and warn user
    if haveIPP && haveIOP
        warning(['Please check whether slices are flipped: ' s.NiftiName]);
    else
        warning(['No orientation/location information found for ' s.NiftiName]);
    end
elseif sign(pos-R(iSL,4)) ~= signSL % same direction?
    R(:,3) = -R(:,3);
end

%% Subfunction: return true if keyword is in s.ImageType
function tf = isType(s, keyword)
typ = tryGetField(s, 'ImageType', '');
tf = ~isempty(strfind(typ, keyword)); %#ok<*STREMP>
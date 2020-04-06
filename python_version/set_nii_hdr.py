import numpy as np
import nibabel as nib

def bitget(number, pos):
    return (number >> pos) & 1

def isfield(s, field):
    """
    Check if string "field" correspond to a key in structure s
    Parameters
    ----------
    s : structure where we search "field" key
    field : key

    Returns
    bool : True if field in s, False otherwise
    -------
    """
    if field in s:
        return True
    else:
        return False

def tryGetField(s, field, dftVal=None):
    """
    Get field if exist, return default value otherwise
    Parameters
    ----------
    s : header in which we want to get the value of the field "field"
    field :  field string value we want to get in s
    dftVal : default return value if field doesn't exist in s

    Returns
    val : value of field "field" in s
    -------
    """
    if field in s:
        val = s[field]
    elif dftVal is not None:
        val = dftVal
    else:
        val = None
    return val

# Subfunction: return true if keyword is in s.ImageType
def isType(s, keyword):
    typ = tryGetField(s, 'ImageType', '')
    tf = not (typ.find(keyword) == -1)  # ok<*STREMP>
    return tf


## Subfunction, Convert 3x3 direction cosine matrix to quaternion
# Simplied from Quaternions by Przemyslaw Baranski
# Retrun quaternion abcd from normalized matrix R (3x3)
def dcm2quat(R):
    proper = np.sign(np.linalg.det(R))
    if proper < 0:
        R[:, 2] = -R[:, 2]
    q = np.sqrt(np.dot(np.array([[1, 1, 1], [1, - 1, - 1], [-1, 1, - 1], [-1, - 1, 1]]), np.diag(R)) + 1) / 2
    if not np.isreal(q[0]):
        q[0] = 0
    mx = np.amax(q)
    ind = np.where(q == mx)
    ind=ind[0][0]
    mx = mx * 4
    if ind[0] == 0:
        q[1] = (R[2, 1] - R[1, 2]) / mx
        q[2] = (R[0, 2] - R[2, 0]) / mx
        q[3] = (R[1, 0] - R[0, 1]) / mx
    elif ind[0] == 1:
        q[0] = (R[2, 1] - R[1, 2]) / mx
        q[2] = (R[0, 1] + R[1, 0]) / mx
        q[3] = (R[2, 0] + R[0, 2]) / mx
    elif ind[0] == 2:
        q[0] = (R[0, 2] - R[2, 0]) / mx
        q[1] = (R[0, 1] + R[1, 0]) / mx
        q[3] = (R[1, 2] + R[2, 1]) / mx
    elif ind[0] == 3:
        q[0] = (R[1, 0] - R[0, 1]) / mx
        q[1] = (R[2, 0] + R[0, 2]) / mx
        q[2] = (R[1, 2] + R[2, 1]) / mx
    if q[0] < 0:
        q = -q  # as MRICron
    return q, proper

def null(a, rtol=1e-5):
    """
   null(A) is an orthonormal basis for the null space of A obtained
   from the singular value decomposition.  That is,  A*Z has negligible
   elements, size(Z,2) is the nullity of A, and Z'*Z = I.
    """
    u, s, v = np.linalg.svd(a)
    rank = (s > rtol * s[0]).sum()
    return rank, v[rank:].T.copy()

def phaseDirection(s):
    phPos = []
    iPhase = []
    fld = 'InPlanePhaseEncodingDirection'
    if isfield(s, fld):
        if s[fld] == 'COL': iPhase = 1  # based on dicm_img(s,0)
        elif s[fld] == 'ROW': iPhase = 0
        else : print(['Unknown ' + fld + ' for ' + s["NiftiName"] + ': ' + s[fld]])
    #TODO : wasn't used in our case (after this, just return the function (case "other manufacturer of the below code)
    # so needs to be checked
    # and needs to implement "csa_header", "bitget"
    # # SIEMENS
    # if isfield(s, 'CSAImageHeaderInfo'):
    #     phPos = False#csa_header(s, 'PhaseEncodingDirectionPositive')  # image ref
    # # GE
    # elif isfield(s, 'ProtocolDataBlock'):
    #     try:  # VIEWORDER "1" == bottom_up
    #         phPos = s["ProtocolDataBlock"]["VIEWORDER"] == '1'
    #     except:
    #         pass
    # # GE
    # elif isfield(s, 'UserDefineData'):
    #     # https://github.com/rordenlab/dcm2niix/issues/163
    #     try:
    #         b = s["UserDefineData"]
    #         i = np.uint16(b[24:26])  # hdr_offset
    #         v = np.single(b[i:i + 4])  # 5.0 to 40.0
    #         if v >= 25.002:
    #             i = i + 76
    #             flag2_off = 777
    #         else:
    #             flag2_off = 917
    #         sliceOrderFlag = bitget(b[i + flag2_off - 1], 2)
    #         phasePolarFlag = bitget(b[i + 49 - 1], 3)
    #         phPos=not (phasePolarFlag ^ sliceOrderFlag)
    #     except:
    #         pass
    # else:
    #     d=''
    #     # Philips
    #     if isfield(s, 'Stack'):
    #         d = s["Stack"]["Item_1"]["MRStackPreparationDirection"][0]
    #         return [phPos, iPhase]
    #     # UIH
    #     elif isfield(s, 'PEDirectionDisplayed'):
    #         d = s["PEDirectionDisplayed"][0]
    #         return [phPos, iPhase]
    #     # Bruker
    #     elif isfield(s, 'Private_0177_1100'):
    #         expr = '(?<=\<\+?)[LRAPSI]{1}(?=;\s*phase\>)'  # <+P;phase> or <P;phase>
    #         d = re.match(expr, str(s.Private_0177_1100), 'once')
    #         id = re.match('LRAPSI', d)
    #         id = id + mod(id,2)*2-1
    #         str = 'LRAPFH'; d = str(id);
    #     # unknown Manufacturer
    #     else:
    #         return [phPos, iPhase]
    #     try:
    #         R = np.reshape(s["ImageOrientationPatient"], (3, 2), 'F')
    #     except:
    #         return [phPos, iPhase]
    #     ixy = np.argmax(np.abs(R), axis=0)
    #     if len(iPhase) == 0:
    #         iPhase = 'RLAPFH'.find(d)
    #         iPhase = np.ceil(iPhase / 2)  # 0/1/2 for RL/AP/FH
    #         iPhase = np.where(ixy == iPhase)  # now 0 or 1
    #     if np.any(d == 'LPH'):
    #         phPos=False  # in dicom ref
    #     elif np.any(d == 'RAF'):
    #         phPos=True
    #     if R[ixy[iPhase], iPhase] < 0: phPos = not phPos  # tricky! in image ref
    #     return [phPos, iPhase]
    return [phPos, iPhase]



## Subfunction: get dicom xform matrix and related info
def xform_mat_manon(s, dim):
    haveIOP = isfield(s, 'ImageOrientationPatient')
    if haveIOP : R = np.reshape(s["ImageOrientationPatient"], (3, 2), order='F')
    else: R = np.array([[1, 0, 0], [0, 1, 0]]).T

    R = np.c_[R, np.cross(R[:, 0], R[:, 1])]  # right handed, but sign may be wrong
    foo = np.abs(R)
    ixyz = np.argmax(foo, axis=1)
    if ixyz[1] == ixyz[0]:
        foo[ixyz[1], 1] = 0
        ixyz[2] = np.argmax(foo[:, 1])
    if np.any(ixyz[2] == ixyz[0:2]): ixyz[2] = np.setdiff1d(np.arange(0, 3), ixyz[0:2])
    iSL = ixyz[2]
    signSL = np.sign(R[iSL, 2])

    try:
        pixdim = s["PixelSpacing"][([1, 0])]
        xyz_unit = 2  # mm
    except:
        pixdim = np.array([1, 1])  # fake
        xyz_unit = 0  # no unit information

    thk = tryGetField(s, 'SpacingBetweenSlices')
    if thk == None: thk = tryGetField(s, 'SliceThickness', pixdim(1))
    pixdim = np.append(pixdim, thk) #warning -> row vector instead of col
    haveIPP = isfield(s, 'ImagePositionPatient')
    if haveIPP : ipp = s["ImagePositionPatient"]
    else : ipp = -(dim.T * pixdim) / 2

    # Next is almost dicom xform matrix, except mosaic trans and unsure slice_dir
    R = np.append(np.dot(R, np.diag(pixdim)), ipp, 1)
    if dim[2] < 2: return ixyz, R, pixdim, xyz_unit

    if isfield(s, 'LastFile') and isfield(s["LastFile"], 'ImagePositionPatient'):
        R[:, 2] = ((s["LastFile"]["ImagePositionPatient"] - R[:, 3].reshape(3, 1)) / (dim[2])).reshape(1, 3)
        thk = np.linalg.norm(R[:, 2])
        if np.abs(pixdim[2] - thk) / thk > 0.01:
            pixdim[2] = thk
            print('slice thickness does not correspond to R')
        return ixyz, R, pixdim, xyz_unit

    #TODO : wasn't used in our case, so needs to be checked and
    # function "csa_header", "asc_header" needs to be implemented
    # elif s["Columns"] > dim[1] and not 'UIH' in s["Manufacturer"][:3]:
    #     R[:, 3] = np.dot(R, np.append(np.ceil(np.sqrt(dim[2]) - 1) * dim[0:2] / 2, [0, 1]))
    #     vec = csa_header(s, 'SliceNormalVector') #don't knox this function called "csa_header"
    #     if not len(vec)==0 # exist for all tested data
    #         if np.sign(vec[iSL]) != signSL :
    #             R[:,2] = -R[:,2]
    #         return ixyz, R, pixdim, xyz_unit
    # if isfield(s, 'CSASeriesHeaderInfo'):  # Siemens both mosaic and regular
    #     print('isfield(s, CSASeriesHeaderInfo')
    #     ori = ['Sag', 'Cor', 'Tra']
    #     ori = ori[iSL]
    #     sNormal = asc_header(s, ['sSliceArray.asSlice[0].sNormal.d' ori])
    #     if asc_header(s, ['sSliceArray.ucImageNumb' ori]):
    #         sNormal = -sNormal
    #     if np.sign(sNormal) != signSL:
    #         R[:,2] = -R[:,2]
    #     if not len(sNormal)==0 :
    #         return ixyz, R, pixdim, xyz_unit
    # pos = []  # volume center we try to retrieve
    # if isfield(s, 'LastScanLoc') and isfield(s, 'FirstScanLocation'):  # GE
    #     pos = (s["LastScanLoc"] + s["FirstScanLocation"]) / 2  # mid-slice center
    #     if iSL < 2: pos = -pos  # RAS convention!
    #     pos = pos - np.dot(R[iSL, 0:2], (dim[0:2].T - 1)) / 2  # mid-slice location
    # if len(pos) == 0 and isfield(s, 'Stack'):  # Philips
    #     ori = ['RL', 'AP', 'FH']
    #     ori = ori[iSL]
    #     pos = tryGetField(s["Stack"]["Item_1"], 'MRStackOffcentre' + ori)
    #     pos = pos - np.dot(R[iSL, 0:2], dim[0:2].T) / 2  # mid-slice location
    # if len(pos) == 0:  # keep right-handed, and warn user
    #     if haveIPP and haveIOP : print('Please check whether slices are flipped: ' + s.NiftiName)
    #     else : print('No orientation/location information found for ' + s.NiftiName)
    # elif np.sign(pos - R[iSL, 3]) != signSL:  # same direction?
    #     R[:, 2] = -R[:, 2]


def normc(M):
    vn = np.sqrt(np.sum(M ** 2, axis=1))
    vn = np.where(vn == 0, 1, vn)
    v = np.true_divide(M, vn)
    return v



def set_nii_hdr(nii, h, pf):
    """
    Complete a blank nifti header (nii file), with computed paramaters (h, pf)

    Parameters
    ----------
    nii : nifti file containing the image and the nifti header who needs to be completed with new computed values
    h : structure of extracted and computed parameters from ismrmrd header
    pf : parameter needed for nifti header completion

    Returns
    nii : completed nifti file with computed parameters
    h : parameters (should have change)
    -------
    """
    #get the nifit file header
    hdr = nii.header
    #get the nifti file image
    nii_empty_img = nii.get_data()
    #get the dimension of the image
    dim = hdr["dim"][1:4]
    #get the number of volumes
    nVol = hdr["dim"][5]

    #Verifies if NumberOfTemporalPositions exists in h and complete it if not
    fld = 'NumberOfTemporalPositions'
    if not fld in h[0] and nVol > 1: h[0][fld] = nVol

    #Compute the transformation matrix: most important feature for nii
    [ixyz, R, pixdim, xyz_unit] = xform_mat_manon(h[0], dim)


    R[0:2, :] = -R[0:2, :]  # dicom LPS to nifti RAS, xform matrix before reorient

    ## Store CardiacTriggerDelayTime
    #TODO : wasn't used in our case, so needs to be checked and function "dicm_hdr" needs to be implemented
    # fld = 'CardiacTriggerDelayTime'
    # if not isfield(h[0], 'CardiacTriggerDelayTimes') and nVol > 1 and isfield(h[0], fld):
    #     if len(h) == 1:
    #         iFrames = np.arange(0, dim[2] * nVol, dim[2])
    #         if isfield(h[0], 'SortFrames') : iFrames = h[0]["SortFrames"][iFrames]
    #         arr = np.empty(0, nVol)
    #         arr[:] = np.nan
    #         s2 = {fld: arr}
    #         s2 = dicm_hdr(h[0], s2, iFrames)
    #         tt = s2[fld]
    #     else:
    #         tt = np.zeros((0, nVol))
    #         inc = int(len(h) / nVol)
    #         for j in range(0, nVol):
    #             tt[j] = tryGetField(h[(j - 1) * inc + 1], fld, 0)
    #     if not (np.all(np.diff(tt) == 0)): h[0]["CardiacTriggerDelayTimes"] = tt

    ## Get EchoTime for each vol
    #TODO : wasn't used in our case, so needs to be checked and function "dicm_hdr" needs to be implemented
    # if not isfield(h[0], 'EchoTimes') and nVol > 1 and isfield(h[0], 'EchoTime'):
    #     if len(h) == 1:
    #         iFrames = np.arange(0, dim[2] * nVol, dim[2])
    #         if isfield(h[0], 'SortFrames'):
    #             iFrames = h[0]["SortFrames"][iFrames]
    #         arr = np.empty(0, nVol)
    #         arr[:] = np.nan
    #         s2 = {'EffectiveEchoTime': arr}
    #         s2 = dicm_hdr(h[0], s2, iFrames)
    #         ETs = s2["EffectiveEchoTime"]
    #     else:
    #         ETs = np.zeros((0, nVol))
    #         inc = len(h) / nVol
    #         for j in range(0, nVol):
    #             ETs[j] = tryGetField(h[(j - 1) * inc + 1], 'EchoTime', 0)
    #     if not (np.all(np.diff(ETs) == 0)):
    #         h[0]["EchoTimes"] = ETs

    #Sets xyz_unit
    hdr["xyzt_units"] = xyz_unit + hdr["xyzt_units"]  # normally: mm (2) + sec (8)

    #Work only on the first header
    s = h[0]
    # Store motion parameters for MoCo series
    #TODO : wasn't used in our case, so needs to be checked
    # if isfield(s, 'RBMoCoTrans') and isfield(s, 'RBMoCoRot') and nVol > 1:
    #     inc = len(h) / nVol
    #     s["RBMoCoTrans"] = np.zeros(nVol, 3)
    #     s["RBMoCoRot"] = np.zeros(nVol, 3)
    #     for j in range(0, nVol):
    #         s["RBMoCoTrans"][j, :] = tryGetField(h[(j - 1) * inc + 1], 'RBMoCoTrans', [0, 0, 0])
    #         s["RBMoCoRot"][j, :] = tryGetField(h[(j - 1) * inc + 1], 'RBMoCoRot', [0, 0, 0])

    # Store FrameReferenceTime: seen in Philips PET
    #TODO : wasn't used in our case, so needs to be checked and needs to create the function "dicm_dict"
    # if isfield(s, 'FrameReferenceTime') and nVol > 1:
    #     inc = len(h) / nVol
    #     vTime = np.zeros(0, nVol)
    #     dict = dicm_dict('', 'FrameReferenceTime')
    #     for j in range (0,nVol):
    #       s2 = dicm_hdr(h{(j-1)*inc+1}.Filename, dict)
    #       vTime[j] = tryGetField(s2, 'FrameReferenceTime', 0)
    #     if vTime(1) > vTime(end) # could also re-read sorted h{i}{1}
    #         vTime = np.flip(vTime)
    #         nii_empty.img = np.flip(nii_empty.img, 4);
    #     s.VolumeTiming = vTime / 1000 # ms to seconds


    # dim_info byte: freq_dim, phase_dim, slice_dim low to high, each 2 bits
    [phPos, iPhase] = phaseDirection(s)  # phPos relative to image in FSL feat!
    if iPhase == 1:
        fps_bits = np.array([1, 4, 16])
    elif iPhase == 0:
        fps_bits = np.array([4, 1, 16])
    else:
        fps_bits = np.array([0, 0, 16])

    ## Reorient if MRAcquisitionType==3D || isDTI && nSL>1
    # If FSL etc can read dim_info for STC, we can always reorient.
    perm = np.argsort(ixyz)  # may permute 3 dimensions in this order
    #Vector used for computing R
    shift_vector = np.array([1, 1, 0.5])

    if ((tryGetField(s, 'MRAcquisitionType', '') == '3D') or s["isDTI"]) and \
            dim[2] > 1 and not (perm == np.arange(0, 3)).all():  # skip if already XYZ order
        R[:, 0:3] = R[:, perm]  # xform matrix after perm
        fps_bits = fps_bits[perm]
        ixyz = ixyz[perm]  # 1:3 after perm
        dim = dim[perm]
        shift_vector = shift_vector[perm]
        pixdim = pixdim[perm]
        hdr["dim"][1:4] = dim
        nii_empty_img = np.transpose(nii_empty_img, perm)
        if isfield(s, 'bvec'):
            s["bvec"] = s["bvec"][:, perm]

    iSL = np.where(fps_bits == 16)
    iPhase = np.where(fps_bits == 4)  # axis index for phase_dim in re-oriented img

    hdr["dim_info"] = np.dot(np.arange(1, 4), np.transpose(fps_bits))  # useful for EPI only
    hdr["pixdim"][1:4] = pixdim  # voxel zize
    hdr["pixdim"][4:8] = [0,0,0,0]
    flp = np.where(R.T.flat[ixyz + [0, 3, 6]] < 0, 1, 0)  # flip an axis if true
    d = np.linalg.det(R[:, 0:3]) * np.prod(1 - flp * 2)  # det after all 3 axis positive
    if (d > 0 and pf["lefthand"]) or (d < 0 and not pf["lefthand"]):
        flp[0] = not flp[0]  # left or right storage

    rotM = np.diag(np.concatenate(((1 - flp * 2), np.array([1]))))  # 1 or -1 on diagnal

    rotM[0:3, 3] = dim * flp  # 0 or dim (initially it was (dim -1) :
    # correlated with the 427th line of this code where I
    # removed the "-1". I think it's because we need to substract 1 on x and y
    # dimensions and 0.5 on z

    R = np.linalg.solve(rotM.T, R.T).T  # xform matrix after flip

    R[:, 3] = R[:, 3] + np.dot(R[:, 0:3], np.transpose(shift_vector))  # This line adjusts values in terms of
    # offsets in x,y and z direction. I saw on many documentations and forums that it's normal to substract [1 1 0.5]*dimpix
    # I havn't verified it for all the possible configurations
    # (positives and negatives rotation angles, etc..) so i'm not really sure about my code,
    # but it works so be suspicious if you encounter issues on R(:,4) and don't hesitate to share your feelings about it !

    for k in range(0, 3):
        if flp[k]: nii_empty_img = np.flip(nii_empty_img, axis=k)
    if flp[iPhase] and phPos:
        phPos = not phPos
    if isfield(s, 'bvec'): s["bvec"][:, flp] = -s["bvec"][:, flp]
    if flp[iSL] and isfield(s, 'SliceTiming'):  # slices flipped
        s["SliceTiming"] = np.flip(s["SliceTiming"])
        sc = hdr["slice_code"]
        if sc > 0:
            hdr["slice_code"] = sc + np.mod(sc, 2) * 2 - 1  # 1<->2, 3<->4, 5<->6

    # sform
    frmCode = np.all([j in s for j in ['ImageOrientationPatient', 'ImagePositionPatient']])
    frmCode = tryGetField(s, 'TemplateSpace', frmCode)
    hdr["sform_code"] = frmCode  # 1: SCANNER_ANAT
    hdr["srow_x"] = R[0, :]
    hdr["srow_y"] = R[1, :]
    hdr["srow_z"] = R[2, :]
    R0 = normc(R[:, 0:3])
    r, sNorm = null(R0[:, (np.setdiff1d(np.arange(0, 3), iSL, assume_unique=True)).T].T)
    if np.sign(sNorm[ixyz[iSL]]) != np.sign(R[ixyz[iSL], iSL]):
        sNorm = -sNorm
 #   shear = np.where((np.linalg.norm(R0[:, iSL[0]] - sNorm) > 0.01), 1, 0)
    R0[:, iSL[0]] = sNorm

    # compute quaternion (qform)
    hdr["qform_code"] = frmCode
    hdr["qoffset_x"] = R[0, 3]
    hdr["qoffset_y"] = R[1, 3]
    hdr["qoffset_z"] = R[2, 3]
    [q, hdr["pixdim"][0]] = dcm2quat(R0)  # 3x3 dir cos matrix to quaternion
    hdr["quatern_b"] = q[1]
    hdr["quatern_c"] = q[2]
    hdr["quatern_d"] = q[3]

    # if shear:
    #     hdr["hdrTilt"] = hdr  # copy all hdr for tilt version
    #     hdr["qform_code"] = 0  # disable qform
    #     gantry = tryGetField(s, 'GantryDetectorTilt', 0)
    #     hdr["hdrTilt"]["pixdim"][iSL[0][0] + 1] = np.linalg.norm(R[0:3, iSL]) * np.cos(gantry * np.pi / 180)
    #     R[0:3, iSL[0][0]] = sNorm * hdr["hdrTilt"]["pixdim"][iSL[0][0] + 1]
    #     hdr["hdrTilt"]["srow_x"] = R[0, :]
    #     hdr["hdrTilt"]["srow_y"] = R[1, :]
    #     hdr["hdrTilt"]["srow_z"] = R[2, :]

    #store some possibly useful info in descrip and other text fields
    str = tryGetField(s, 'ImageComments', '')
    if isType(s, '\MOCO\\') : str = ''  # seless for MoCo
    foo = tryGetField(s, 'StudyComments')
    if foo is not None : str = str + ';' + foo
    str = str + ';' + s["Manufacturer"].split(' ', 1)[0]
    foo = tryGetField(s, 'ProtocolName')
    if foo is not None : str = str + ";" + foo
    hdr["aux_file"] = str  # char[24], info only

    if not iPhase == []:
        if phPos == []:
            pm = '?'
            b67 = 0
        elif phPos:
            pm = ''
            b67 = 1
        else:
            pm = '-'
            b67 = 2
        hdr["dim_info"] = hdr["dim_info"] + b67 * 64
        axes = 'xyz'  # actually ijk
        phDir = pm + axes[iPhase[0][0]]
        s["UnwarpDirection"] = phDir

    #Slope and intercept: apply to img if no rounding error
    #TODO :  wasn't used in our case, so needs to be checked
    # sclApplied = tryGetField(s, 'ApplyRescale', False)
    # if isfield(s, 'RescaleSlope') and isfield(s, 'RescaleIntercept') and not sclApplied:
    #     slope = tryGetField(s, 'RescaleSlope', 1)
    #     inter = tryGetField(s, 'RescaleIntercept', 0)
    #     if isfield(s, 'MRScaleSlope'):  # Philips: see PAR file for detail
    #         inter = inter / (slope * float(s["MRScaleSlope"]))
    #         slope = 1 / float(s["MRScaleSlope"])
    #     val = list(np.array([np.max(nii.img[:]), np.min(nii.img)[:]])* slope + inter)
    #     val.sort()
    #     dClass =class(nii.img)
    #         if isa(nii.img, 'float') or (
    #                 mod(slope, 1) == 0 and mod(inter, 1) == 0 and val(1) >= intmin(dClass) and val(2) <= intmax(
    #                 dClass)):
    #             nii.img = nii.img * slope + inter
    #         else:
    #             hdr["scl_slope"] = slope
    #             hdr["scl_inter"] = inter
    # elif sclApplied and isfield(s, 'MRScaleSlope'):
    #     slope = tryGetField(s, 'RescaleSlope', 1) * s["MRScaleSlope"]
    #     nii.img = nii.img / slope

    h[0] = s

    #return the new nifti file
    nii = nib.Nifti1Image(nii_empty_img, None, header=hdr)

    return [nii, h]

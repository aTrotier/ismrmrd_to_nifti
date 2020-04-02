import sys
sys.path.append('/mnt/Lourd/Work/Code/Valery/Code_Matlab/ismrmrd_to_nifti')

from python_version import set_nii_hdr as tools
from python_version import load_ismrmrd_ifft3d_reconstruction as h5reco
from python_version import extract_ismrmrd_parameters_from_headers as param
from python_version import flip_image as fi

import numpy as np
import nibabel as nib
import os
import subprocess
import re
from tkinter import filedialog
from tkinter import *


############################# SELECT A CONVERSION TYPE ###############################

## NEW ISMRMRD VERSION
# 0  :  IF conversion already done with the new ismrmrd version

## OLD ISMRMRD VERSION
# 1 :   (TODO) IF conversion from .dat to .h5 needed (Export from : dat -> ismrmrd -> nifti)
#       you should use the dedicated parameter maps for the siemens_to_ismrmrd
#       conversion in order to have access to all the parameters used for nifti
#       reconstruction

# 2 :   IF conversion already done with our parameter maps : Select .h5 file 
#       (.dat is not required if you do the conversion with our parameter cards too)

# 3  :  IF conversion already done but not with our parameter maps (need .H5 and .dat)
 

conversion_type = 3


#################### SELECT A OUPUT PATH ###############################
# If output_path = '' : default path = data_path/filename+_ismrmrd_to_nifti_version_python.nii

output_path = ''













 
#####################################################################################
############################## DON'T CHANGE ANYTHING FROM HERE ######################

## CASE 0 :
if conversion_type == 0 :
    root = Tk()
    root.filename =  filedialog.askopenfilename(initialdir = "/",title = "Select ISMRMRD file (*.h5)",filetypes = [("ISMRMRD files","*.h5")])
    print (root.filename)
    filename  = root.filename

## CASE 1 : TODO
if conversion_type == 1 :
    xmlFile = '/home/.../ismrmrd_to_nifti/parameterMaps/IsmrmrdParameterMap_Siemens_Table.xml'
    xlsFile = '/home/.../ismrmrd_to_nifti/parameterMaps/IsmrmrdParameterMap_Siemens_Table.xsl'
    
    filename_raw = '.dat'

## CASE 2 :
if conversion_type == 2 :
    root = Tk()
    root.filename =  filedialog.askopenfilename(initialdir = "/",title = "Select ISMRMRD file (*.h5)",filetypes = [("ISMRMRD files","*.h5")])
    filename = root.filename

## CASE 3 : 
if conversion_type == 3 :
    root = Tk()
    root.filename =  filedialog.askopenfilename(initialdir = "/",title = "Select ISMRMRD file (*.h5)",filetypes = [("ISMRMRD files","*.h5")])
    filename = root.filename
    root.filename =  filedialog.askopenfilename(initialdir = "/",title = "Select RAW SIEMENS file (*.dat)",filetypes = [("RAW SIEMENS files","*.dat")])
    filename_raw = root.filename


## Load selected file
[pathstr, name] = os.path.split(filename)
[name, ext]=os.path.splitext(name)

## Reconstruct image (z,y,x)
[head, hdr, img_scaled] = h5reco.load_ismrmrd_ifft3d_reconstruction(filename)

## Correct "position" fields with the table position offset parameter
if conversion_type == 1 or conversion_type == 2:
    #Find indices of GlobalTableProsTra fields
    idx_TablePosTra = [i for i in range(len(hdr.userParameters.userParameterLong)) if  hdr.userParameters.userParameterLong[i].name == 'GlobalTablePosTra']
    if idx_TablePosTra == []:
        print('Missing GlobalTablePosTra parameters in the ismrmrd dataset. Check if you use the parameter maps in the folder')
        print("Exit with code 0.")
        sys.exit(0)
    else :
        table_position_offset = int(hdr.userParameters.userParameterLong[idx_TablePosTra[0]].value_)
elif conversion_type == 3 :
    # Get table_position_offset : really important ! get value of the "0" position of the table in our case (SIEMENS AREA)..
    # to complete and confirm for other systems
    # need GlobalTablePosTra (table position offset) contained by siemens raw data
    GlobalTablePosTra =  subprocess.check_output('grep -a "GlobalTablePosTra" '+ filename_raw, shell=True).decode()
    GlobalTablePosTra=re.findall(' \-*\d+', GlobalTablePosTra)
    table_position_offset = int(GlobalTablePosTra[0])
else :
    table_position_offset = 0 
head.position[2] += table_position_offset

## Create parameters for set_nii_hdr et xform_mat
h = param.extract_ismrmrd_parameters_from_headers(head, hdr)

## Get crop image, flip and rotationate to match with true Nifti image
img = fi.flip_image(img_scaled[0])

## Create nii struct based on img
nii_empty = nib.Nifti1Image(img, np.eye(4))
pf = {"lefthand":0}
h2 =[]
h2.append(h)

## Compute nifti parameters
hdr=nii_empty.header
[nii_filled, h3] = tools.set_nii_hdr(nii_empty, h2, pf)

## Save image in nifti format
if output_path == '' :
    output_path = os.path.join(pathstr +'/'+ name +'_ismrmrd_to_nifti_version_python.nii')
nib.save(nii_filled, output_path)
print("Nifti well saved")



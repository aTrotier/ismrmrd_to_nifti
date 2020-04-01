from ismrmrd_to_nifti import set_nii_hdr as tools
from ismrmrd_to_nifti import load_ismrmrd_ifft3d_reconstruction as h5reco
from ismrmrd_to_nifti import extract_ismrmrd_parameters_from_headers as param
from ismrmrd_to_nifti import flip_image as fi
import numpy as np
import nibabel as nib
import os



##IF conversion from .dat to .h5 needed (Export from : dat -> ismrmrd -> nifti)

#you should use the dedicated parameter maps for the siemens_to_ismrmrd
#conversion in order to have access to all the parameters used for nifti
#reconstruction :
xmlFile = '/home/.../ismrmrd_to_nifti/parameterMaps/IsmrmrdParameterMap_Siemens_Table.xml'
xlsFile = '/home/.../ismrmrd_to_nifti/parameterMaps/IsmrmrdParameterMap_Siemens_Table.xsl'


#ELSE : Select .h5 file (.dat is not required if you do the conversion with our parameter cards too)
#filename = '/home/manondesclides/Data/Dataset5/FID/meas_MID00033_FID13782_gre3D_2_2_tranversal_A2P.h5'
filename  = '/media/manon/Lourd/Work/Data/Dataset7/Sheep_Ablation_2/10/meas_MID00136_FID16393_CS_ANGIO_1_3ISO.h5'


# Load selected file
[pathstr, name] = os.path.split(filename)
[name, ext]=os.path.splitext(name)
#reconstruct image (z,y,x)
[head, hdr, img_scaled] = h5reco.load_ismrmrd_ifft3d_reconstruction(filename)

# Get table_position_offset : really important ! get value of the "0" position of the table in our case (SIEMENS AREA)..
# to complete and confirm for other systems
# need GlobalTablePosTra (table position offset) contained by siemens raw data
#GlobalTablePosTra =  subprocess.check_output('grep -a "GlobalTablePosTra" '+ filename_raw, shell=True).decode()
#GlobalTablePosTra=re.findall(' \-*\d+', GlobalTablePosTra)
#table_position_offset=int(GlobalTablePosTra[0])

# Create parameters for set_nii_hdr et xform_mat
h = param.extract_ismrmrd_parameters_from_headers(head, hdr)

#Get crop image, flip and rotationate to match with true Nifti image
img = fi.flip_image(img_scaled[0])

#Convert and save our new image in Nifti format
#create nii struct based on img
nii_empty = nib.Nifti1Image(img, np.eye(4))
pf = {"lefthand":0}
h2 =[]
h2.append(h)

#Save image in nifti format
hdr=nii_empty.header
[nii_filled, h3] = tools.set_nii_hdr(nii_empty, h2, pf)

# #save the image
nib.save(nii_filled, os.path.join(pathstr +'/'+ name +'_ismrmrd_to_nifti_version_python.nii'))
print("Nifti well saved")



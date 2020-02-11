# ismrmrd_to_nifti

New beta version of ismrmrd_to_nifti ("manon" folder)

## Description 

This new version includes :
- a **'table_position_offset'** parameter which is important to convert well ismrmrd to nifti images (in case of SIEMENS images). This parameter was already in ismrmrd header as 'table_position_offset' but the value wasn't usable beacause it was not the absolute offset, but just the position of the table (need a pair of values : origin AND position to compute the offset). The origin value is stocked in Siemens RAW data and the absolute difference between origin and position too. We get this last value by making a grep in the raw data file associated with the '.h5' to convert. 
This offset is added to the 'z' patient position coordinate.

- little changes in the conversion formula : we hilighted issues on transformation matrix concerning the position of the image. Those little differences between the true Nifti image and the converted one were due to a bad way to deal with matrices dimensions. Indeed, it's necessary to substract 1 voxel dimension on x and y coordinates AND 1/2 slice dimension on z coordinate to center the image. 
This part is a bit more dark and real explanations are missing, so if you know why it's necessary, don't hesitate to share your ideas and change the formula !



Thoses changes have been tested on SIEMENS AREA system, with Sagittal, Coronal and Transversal views, for all read directions and rotations (positives and negatives but not for all cases) in HFS.

If you need our dataset to test the code, tell us !


## Usage

This beta version is available in the "manon" folder.
Her you will find :

- *"main_conversion.m"* : main code to **run conversion** from ".dat"/".h5" to Nifti ".nii"
- *"load_ismrmrd_ifft3d_reconstruction.m"* : function to load ".h5" file and reconstruct images
- *"extract_ismrmrd_parameters_from_headers.m"* : function to compute and stock parameters from ismrmrd format needed by *"set_nii_hdr.m"*
- *"set_nii_hdr.m"* : function to compute and create nifti structure from computed parameters in *"extract_ismrmrd_parameters_from_headers.m" *
- *"flip_image.m"* : function to flip the reconstructed image and fit with the true Nifti one.
- *"comparison_test.m"* : main code to **compare two nifti files** (usefull to compare true and converted nifti files)


### To convert data : 
1) Open **"main_conversion.m"** :

1) In "%% Add path to needed librairies", change paths for your computer

```matlab
addpath('/usr/local/share/ismrmrd/matlab/')
addpath('/home/path_to_ismrmrd_to_nifti')
addpath('/home/path_to_ismrmrd_to_nifti/xiangruili-dicm2nii-b76a158')
addpath('/home/path_to_ismrmrd_to_nifti/manon')

```
3)
    ---> If you need to **convert your file to Nifti from .dat**, uncomment the part "%% IF conversion from .dat to .h5 needed (Export of MP2RAGE-CS_WIP siemens from : dat -> ismrmrd -> nifti)
" and set the path *filename_raw* to your ".dat" file.

    ---> If you need to **convert your file to Nifti from .h5** uncomment the part "%% ELSE : Select .h5 file (need .dat too)" and set the path *filename* to your ".h5" file and the *filename_raw* to your ".dat" file.

3) Change the last line of the code to save in a different name or path 
```matlab
nii_tool('save', nii2, [pathstr '/' name '_ismrmrd_to_nifti_version.nii'], rst3D);
```
if you want"

4) Run the code

### To compare real Nifti format and the converted one : 
1) Open **"comparison_test.m"**
1) Set the *true_filename* to the true Nifti dataset
2) Set the *converted_filename* to the converted Nifti dataset
3) Run and compare



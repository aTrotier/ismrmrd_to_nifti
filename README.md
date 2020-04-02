# ismrmrd_to_nifti

New beta version of ismrmrd_to_nifti

## Data used for test :
The dataset we worked on is available here : https://zenodo.org/record/3674622#.Xk031Y7Pw5k

## Usage

Her you will find :

- *"main_conversion.m"* : main code to **run conversion** from ".dat"/".h5" to Nifti ".nii"
- *"load_ismrmrd_ifft3d_reconstruction.m"* : function to load ".h5" file and reconstruct images
- *"load_reconstructed_data.m"* : function to load ".mat" file of ismrmrd images
- *"extract_ismrmrd_parameters_from_headers.m"* : function to compute and stock parameters from ismrmrd format needed by *"set_nii_hdr.m"*
- *"set_nii_hdr.m"* : function to compute and create nifti structure from computed parameters in *"extract_ismrmrd_parameters_from_headers.m" *
- *"flip_image.m"* : function to flip the reconstructed image and fit with the true Nifti one.
- *"crop_image.m"* : function to crop the reconstructed image and fit with the reconstructed data size.
- *"comparison_test.m"* : main code to **compare two nifti files** (usefull to compare true and converted nifti files)


### To convert data : 
1) Open **"main_conversion.m"** :

1) In "%% Add path to needed librairies", change paths for your computer

```matlab
addpath('/usr/local/share/ismrmrd/matlab/')
addpath('/home/path_to_ismrmrd_to_nifti')
addpath('/home/path_to_ismrmrd_to_nifti/xiangruili-dicm2nii-b76a158')

```
2) Choose a type of conversion (*conversion_type*) depending on your data type :
    For **NEW VERSION of siemens_to_ismrmrd** (git commit > May 15th, 2019) :
    - **0**  :  If conversion already done with the new siemens_to_ismrmrd version (git commit > May 15th, 2019) : select the .h5 or .mat corresponding to your data file.

    For **OLD VERSION of siemens_to_ismrmrd** (git commit < May 15th, 2019) :
    - **1** :   If conversion from .dat to .h5 needed (Export from : dat -> ismrmrd -> nifti) you should use the dedicated parameter maps for the siemens_to_ismrmrd conversion in order to have access to all the parameters used for nifti reconstruction (parameter maps in *parameterMaps* folder)
    - **2** :   If conversion already done with our parameter maps : Select .h5 file or .mat corresponding to your data file.
    - **3** :   If conversion already done but **NOT** with our parameter maps : Select .h5 file or .mat **AND** .dat corresponding to your data file.

3) Change the "*ouput_path*" to save in a different name or path (by default : *data_path/filename+_ismrmrd_to_nifti_version_matlab.nii* for matlab version and *data_path/filename+_ismrmrd_to_nifti_version_python.nii* for python version).


4) Run the code

### To compare real Nifti format and the converted one : 
1) Open **"comparison_test.m"**
1) Set the *true_filename* to the true Nifti dataset
2) Set the *converted_filename* to the converted Nifti dataset
3) Run and compare


## Description 

This new version includes :
- a Matlab code to convert ISMRMRD data files into Nifti format from ".h5" or ".mat" files (available in *matlab_version* folder);
- a Python code to convert ISMRMRD data files into Nifti format from ".h5" or ".mat" files (available in *python_version* folder);
- this version allows conversion for ismrmrd data converted before the addition of the table_position_offset parameter to the patient position (siemens_to_ismrmrd before May 15th 2019) and after by choosing the data type
- little changes in the conversion formula : we hilighted issues on transformation matrix concerning the position of the image. Those little differences between the true Nifti image and the converted one were due to a bad way to deal with matrices dimensions. Indeed, it's necessary to substract 1 voxel dimension on x and y coordinates AND 1/2 slice dimension on z coordinate to center the image. 
This part is a bit more dark and real explanations are missing, so if you know why it's necessary, don't hesitate to share your ideas and change the formula !



Thoses changes have been tested on SIEMENS AREA system, with Sagittal, Coronal and Transversal views, for all read directions and rotations (positives and negatives but not for all cases) in HFS.

If you need our dataset to test the code, tell us !

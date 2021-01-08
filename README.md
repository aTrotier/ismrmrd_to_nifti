# ismrmrd_to_nifti
## Overview

This repository contains two convertors :
* **img_h5_to_nii** -> reconstructed images in .h5 from gadgetron to nifti
* **gadgetron_to_nii** -> reconstruct nifti from a matlab gadget called in gadgetron (work for bucket)
* **raw_h5_to_nii** -> dumped .h5 rawdata to nifti (and associated function that can be used to implement the export in a matlab/python gadgets)

## img_h5_to_nii

For now only the matlab convertor is implemented (feel free to dev the python one).
The main function is : *ismrmrd_img_to_nifti.m*
## raw_h5_to_nii
Implementation in matlab and python. Check the [README.md](./raw_h5_to_nii/README.md) in the subfolder.



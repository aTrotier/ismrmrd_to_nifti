#!/usr/bin/env python3

import sys
sys.path.append('/opt/src')

from python_version import set_nii_hdr as tools
from python_version import load_ismrmrd_ifft3d_reconstruction as h5reco
from python_version import extract_ismrmrd_parameters_from_headers as param
from python_version import flip_image as fi

import numpy as np
import nibabel as nib
import os
import subprocess
import re
import getopt


def main(argv):
    inputfile = ''
    outputfile = ''
    try:
        opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
    except getopt.GetoptError:
        print('test.py -i <inputfile> -o <outputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('ismrmrd_to_nifti-python -i <inputfile> -o <outputfile>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg


    ## Load selected file
    filename = inputfile
    
    [pathstr, name] = os.path.split(filename)
    [name, ext]=os.path.splitext(name)

    if outputfile == '' :
        output_path = os.path.join(pathstr +'/'+ name +'_ismrmrd_to_nifti_version_python.nii')
    else :
        output_path = outputfile 
    print('Input file : ', filename)
    print('Output file : ', output_path)

    ## Reconstruct image (z,y,x)
    [head, hdr, img_scaled] = h5reco.load_ismrmrd_ifft3d_reconstruction(filename)

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
    [nii_filled, h3] = tools.set_nii_hdr(nii_empty, h2, pf)

    ## Save image in nifti format
    nib.save(nii_filled, output_path)

    print("Nifti well saved")


if __name__ == "__main__":
   main(sys.argv[1:])


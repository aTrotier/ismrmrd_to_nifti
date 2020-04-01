
import numpy as np


def flip_image(img):
    """

    Parameters
    ----------
    img : image to flip and rotate to have the same orientation than the real nifti image

    Returns
    -------
    flipped_img : the flipped image

    """
    #permute all dimensions to have yxz order
    img = img.transpose(1, 2, 0)
    flipped_img = img

    return flipped_img
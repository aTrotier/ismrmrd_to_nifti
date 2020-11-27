import numpy as np
import sys

def extract_ismrmrd_parameters_from_headers(head, hdr):
    """

    Parameters
    ----------
    head : read imsrmrd dataset head (dataset.head)
    hdr : deserialized ismrmrd xml dataset file
    Returns
    -------
    function that returns a structure h which contains parameters needed for conversion function
    """

    #get encoding parameters from hdr file
    enc = hdr.encoding[0]
    # Matrix size
    rNx = enc.reconSpace.matrixSize.x
    rNy = enc.reconSpace.matrixSize.y
    rNz = enc.reconSpace.matrixSize.z
    # Field of View
    rFOVx = enc.reconSpace.fieldOfView_mm.x
    rFOVy = enc.reconSpace.fieldOfView_mm.y
    rFOVz = enc.reconSpace.fieldOfView_mm.z

    # Create structure h and fill it with appropriate values
    h = {"NumberOfTemporalPositions": 1,
         "MRAcquisitionType": "3D",
         "ImageOrientationPatient": np.concatenate((-np.array(head.phase_dir), -np.array(head.read_dir)), axis=0),
         "SpacingBetweenSlices": rFOVz / rNz,
         "SliceThickness": rFOVz / rNz,
         "PixelSpacing": np.array([[rFOVx / rNx], [rFOVy / rNy]])}

    if (hdr.encoding[0].reconSpace.matrixSize.z == 0):
        hdr.encoding[0].reconSpace.matrixSize.z = 1


    #ImagePositionPatient of first slice (right corner position of first and last slice)
    #Attention peut-être que faut changer les +/-
    #en fonction de HFS dans : hdr.measurementInformation.patientPosition
    corner_first_slice = np.zeros((3, 1), dtype=float)
    corner_last_slice = np.zeros((3, 1), dtype=float)

    #Compute the corner position (ImagePositionPatient) of the first slice
    corner_first_slice[0] = head.position[0] \
                            + rFOVx / 2 * head.read_dir[0] \
                            + rFOVy / 2 * head.phase_dir[0] \
                            + rFOVz / 2 * head.slice_dir[0]
    corner_first_slice[1] = head.position[1] \
                            + rFOVx / 2 * head.read_dir[1] \
                            + rFOVy / 2 * head.phase_dir[1] \
                            + rFOVz / 2 * head.slice_dir[1]
    corner_first_slice[2] = head.position[2] \
                            + rFOVx / 2 * head.read_dir[2] \
                            + rFOVy / 2 * head.phase_dir[2] \
                            + rFOVz / 2 * head.slice_dir[2]
    h.update({"ImagePositionPatient": corner_first_slice})

    #Compute the corner position (ImagePositionPatient) of the last slice
    corner_last_slice[0] = head.position[0] \
                           + rFOVx / 2 * head.read_dir[0] \
                           + rFOVy / 2 * head.phase_dir[0] \
                           - rFOVz / 2 * head.slice_dir[0]
    corner_last_slice[1] = head.position[1] \
                           + rFOVx / 2 * head.read_dir[1] \
                           + rFOVy / 2 * head.phase_dir[1] \
                           - rFOVz / 2 * head.slice_dir[1]
    corner_last_slice[2] = head.position[2] \
                           + rFOVx / 2 * head.read_dir[2] \
                           + rFOVy / 2 * head.phase_dir[2] \
                           - rFOVz / 2 * head.slice_dir[2]

    h.update({"LastFile": {"ImagePositionPatient": corner_last_slice}})

    h.update({"Manufacturer": hdr.acquisitionSystemInformation.systemVendor})
    h.update({"Columns": hdr.encoding[0].reconSpace.matrixSize.y})
    h.update({"NiftiName": hdr.measurementInformation.protocolName})
    h.update({"InPlanePhaseEncodingDirection": 'COL'})

    return h
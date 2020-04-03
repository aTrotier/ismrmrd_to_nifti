import os
import ismrmrd
import ismrmrd.xsd
import numpy as np

from ismrmrdtools import transform

def load_ismrmrd_ifft3d_reconstruction(filename):
    """
    Load .h5 file, read header (head) and reconstruct images

    Parameters
    ----------
    filename : path to the .h5 file

    Returns
    -------
    head : read imsrmrd dataset head (dataset.head)
    hdr : deserialized ismrmrd xml dataset file
    img_scaled : reconstructed image

    """

    if not os.path.isfile(filename):
        print("%s is not a valid file" % filename)
        raise SystemExit
    dset = ismrmrd.Dataset(filename, 'dataset', create_if_needed=False)

    #Read some fields from the XML header
    hdr = ismrmrd.xsd.CreateFromDocument(dset.read_xml_header())
    #get encoding and reconstruction information
    enc = hdr.encoding[0]
    # Matrix size
    eNx = enc.encodedSpace.matrixSize.x
    eNy = enc.encodedSpace.matrixSize.y
    eNz = enc.encodedSpace.matrixSize.z
    rNx = enc.reconSpace.matrixSize.x
    rNy = enc.reconSpace.matrixSize.y

    # Number of Slices, Reps, Contrasts, etc.
    #We have to wrap the following in a if/else because a valid xml header may
    #not have an entry for some of the parameters
    ncoils = hdr.acquisitionSystemInformation.receiverChannels
    if enc.encodingLimits.slice != None:
        nslices = enc.encodingLimits.slice.maximum + 1
    else:
        nslices = 1

    if enc.encodingLimits.repetition != None:
        nreps = enc.encodingLimits.repetition.maximum + 1
    else:
        nreps = 1

    if enc.encodingLimits.contrast != None:
        ncontrasts = enc.encodingLimits.contrast.maximum + 1
    else:
        ncontrasts = 1


    # Loop through the acquisitions looking for noise scans
    firstacq = 0
    for acqnum in range(dset.number_of_acquisitions()):
        acq = dset.read_acquisition(acqnum)

        # TODO: Currently ignoring noise scans
        if acq.isFlagSet(ismrmrd.ACQ_IS_NOISE_MEASUREMENT):
            print("Found noise scan at acq ", acqnum)
            continue
        else:
            firstacq = acqnum
            print("Imaging acquisition starts acq ", acqnum)
            break

    # Initialiaze a storage array
    all_data = np.zeros((nreps, ncontrasts, nslices, ncoils, eNz, eNy, rNx), dtype=np.complex64)

    # Loop through the rest of the acquisitions and stuff
    for acqnum in range(firstacq, dset.number_of_acquisitions()):
        acq = dset.read_acquisition(acqnum)
        head = acq.getHead()

        # TODO: this is where we would apply noise pre-whitening

        #padd if acquisition data is not complete (padding)
        if acq.data.shape[1]<eNx :
            x0=int((eNx - acq.data.shape[1]) / 2)
            zeros = np.zeros((acq.data.shape[0], x0))
            padded_acq_data = np.append(np.append(zeros, acq.data, axis=1), zeros, axis=1)
            acq.resize(eNx, acq.active_channels, acq.trajectory_dimensions)
            acq.data[:]=padded_acq_data

        # Remove oversampling if needed
        if eNx != rNx:
            #xline = transform.transform_kspace_to_image(acq.data, [1])
            xline = transform.transform_kspace_to_image(acq.data, dim=(1,), img_shape=(eNx,))
            x0 = int((eNx - rNx) / 2)
            x1 = int((eNx - rNx) / 2 + rNx)
            xline = xline[:, x0:x1]
            acq.resize(rNx, acq.active_channels, acq.trajectory_dimensions)
            acq.center_sample = int(rNx / 2)
            # need to use the [:] notation here to fill the data
            acq.data[:] = transform.transform_image_to_kspace(xline, dim=(1,), k_shape=(rNx,))

        # Stuff into the buffer
        rep = acq.idx.repetition
        contrast = acq.idx.contrast
        slice = acq.idx.slice
        y = acq.idx.kspace_encode_step_1
        z = acq.idx.kspace_encode_step_2
        all_data[rep, contrast, slice, :, z, y, :] = acq.data

    # Reconstruct images
    images = np.zeros((nreps, ncontrasts, nslices, eNz, rNy, rNx), dtype=np.float32)
    img_scaled = []
    for rep in range(nreps):
        for contrast in range(ncontrasts):
            for slice in range(nslices):
                # FFT
                if eNz > 1:
                    # 3D
                    im = transform.transform_kspace_to_image(all_data[rep, contrast, slice, :, :, :, :], [1, 2, 3])
                else:
                    # 2D
                    im = transform.transform_kspace_to_image(all_data[rep, contrast, slice, :, 0, :, :], [2, 3])

                if eNy != rNy:
                    x0 = int((eNy - rNy) / 2)
                    x1 = int((eNy - rNy) / 2 + rNy)
                    im = im[:,:,x0:x1, :]

                # Sum of squares
                im = np.sqrt(np.sum(np.abs(im) ** 2, 0))

                # Stuff into the output
                if eNz > 1:
                    # 3D
                    images[rep, contrast, slice, :, :, :] = im
                else:
                    # 2D
                    images[rep, contrast, slice, 0, :, :] = im

                img_scaled.append(im)

    dset.close()

    return [head, hdr, img_scaled]


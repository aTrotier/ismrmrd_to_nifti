function [head, hdr, img_scaled] = load_reconstructed_data(filename)
    %load data .mat
    b = load(filename)
    %get reconstrcution parameters 
    rec_Nx = b.h.encoding.reconSpace.matrixSize.x;
    rec_Ny = b.h.encoding.reconSpace.matrixSize.y;
    rec_Nz = b.h.encoding.reconSpace.matrixSize.z;

    img_scaled = crop_image({b.img_scale},rec_Nx, rec_Ny, rec_Nz);
    
    hdr=b.h;

    head = b.data.headers;
    
end


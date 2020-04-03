function img_scaled = crop_image(reconImages,rec_Nx, rec_Ny, rec_Nz)


img_scaled = reconImages;
for i=1:size(reconImages)

    if (size(img_scaled{i},1)>rec_Nx)
        shift_x=round((size(img_scaled{i},1)-rec_Nx)/2);
        new_img=img_scaled{i}(shift_x+1:end-shift_x,:,:);
        clear img_scaled{i}
        img_scaled{i}=new_img;
        disp('cropping along x direction due to oversampling')
    end
    if (size(img_scaled{i},2)>rec_Ny)
        shift_y=round((size(img_scaled{i},2)-rec_Ny)/2);
        new_img=img_scaled{i}(:,shift_y+1:end-shift_y,:);
        clear img_scaled
        img_scaled{i}=new_img;
        disp('cropping along y direction due to oversampling')
    end
    if (size(img_scaled{i},3)>rec_Nz)
        shift_z=round((size(img_scaled{i},3)-rec_Nz)/2);
        new_img=img_scaled{i}(:,:,shift_z+1:end-shift_z);
        clear img_scaled
        img_scaled{i}=new_img;
        disp('cropping along z direction due to oversampling')
    end

end


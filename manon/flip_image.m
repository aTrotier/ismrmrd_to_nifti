function flipped_img = flip_image(img)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

%get only first image (here its the only one)
img = permute(img, [2 1 3]);

%OK pour coronal, sagittal et transversal:
img=flip(img,3);
img=flip(img,2);
img = rot90(rot90(img,1),1);
flipped_img = img;

end


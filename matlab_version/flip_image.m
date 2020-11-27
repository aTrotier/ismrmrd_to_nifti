function flipped_img = flip_image(img)

%permute first and second dimension to have yxz order
img = permute(img, [2 1 3 4]);
flipped_img = img;

end


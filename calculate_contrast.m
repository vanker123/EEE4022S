function Image_Contrast = calculate_contrast(IsarImage)

I_temp = abs(IsarImage).^2;
mean_I_temp= mean(I_temp,'all');

IC_Numerator=sqrt(mean((I_temp-mean_I_temp).^2,'all'));
Image_Contrast=IC_Numerator/mean_I_temp;

end 


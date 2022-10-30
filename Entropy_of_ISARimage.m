function Entropy_ISARimage = Entropy_of_ISARimage(ISARimage)


Normalised_power = (abs(ISARimage).^2)/(sum(sum(abs(ISARimage).^2)));
Entropy_ISARimage = -sum(sum(Normalised_power.*log(Normalised_power)));

end 
function Range_Alignment = aligned_range(HRR_profiles)

unaligned_Range= abs(HRR_profiles);
x=unaligned_Range(1,:);   % 1st range profile

len=length(unaligned_Range);
si=size(unaligned_Range);

for i = [1:si(1)]
    unaligned= unaligned_Range(i,:);
    [shift,lg]=xcorr(x,unaligned);
    [max_val max_ind]= max(shift);
    BinShiftsNeeded = max_ind - length(x);
    Shifted_Profile=circshift(unaligned,BinShiftsNeeded);
    AlignedProfiles_Matrix(i, :) = Shifted_Profile;
    BinShiftsNeeded_vector(i) = BinShiftsNeeded;
end
%figure; plot(BinShiftsNeeded_vector)
%xlabel("Range Profile")
%ylabel("Integer bin shifts")
%figure; imagesc(20*log10(AlignedProfiles_Matrix));
%axis xy;
%colormap('jet');
%set(gcf,'color','w')
%set(gca,'FontSize',11)
%ylabel('Profile Number','fontsize',14);
%xlabel('Range (m)','fontsize',14);
%title('HRR profiles','fontsize',14);

 %% Haywood Alignment algorithm

NumberOfProfiles = size(unaligned_Range,1); % Assumption: each row is a range profile
RangeProfilesVector = 1:1:NumberOfProfiles;
order=1;
BestFit=polyfit(RangeProfilesVector,BinShiftsNeeded_vector, order);
BinShiftsNeededNonInteger = polyval(BestFit, RangeProfilesVector);

%figure; 
%plot(BinShiftsNeeded_vector,'-o');
%hold on;
%plot(RangeProfilesVector, BinShiftsNeededNonInteger,'-*r');
%xlabel("Range Profile")
%ylabel("Non integer bin shift")
%set(gcf,'color','w')

HaywoodAligned=zeros(si(1),si(2));

number_of_rows=size(unaligned_Range,1);
number_of_columns=size(unaligned_Range,2);

for j = 1:number_of_rows

    N = length(unaligned_Range);
    n = 0:1:(N-1);
    k = transpose(BinShiftsNeededNonInteger);
    
    Phi_vector = exp(-1i*2*pi*n.*k(j)/N);
    Range_Profile_shifted = ifft(fft(HRR_profiles(j,:)).*Phi_vector);
    HaywoodAligned(j,:)=Range_Profile_shifted;  
    
end 

%figure;
Range_Alignment=HaywoodAligned;
%imagesc(20*log10(abs(HaywoodAligned)));
%xlabel('Range (m)');
%ylabel("doppler frequency")
%axis xy;
%colormap('jet');
%colorbar
%set(gca,'FontSize',14)
%set(gcf,'color','w')
%ylabel('Profile Number','fontsize',14);
%xlabel('Range (m)','fontsize',14);
%title('HRR profiles','fontsize',14);
end 


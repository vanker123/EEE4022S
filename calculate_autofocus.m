function [Plot_IsarImage,IsarImage] = calculate_autofocus(Aligned_HRR_profiles,xAxis,yAxis)

Shifted_profile=abs(Aligned_HRR_profiles);
number_of_rows=size(Aligned_HRR_profiles,1);
number_of_columns=size(Aligned_HRR_profiles,2);

%Amplitude variance to find dominant scatterer 
%computing variance of aligned data

var_rangebins = var(Shifted_profile,1);
%figure;
%plot(var_rangebins)
%xlabel('Range(m)','fontsize',15)
%ylabel('Amplitude variance','fontsize',15)
set(gcf,'color','w')

sum_profile= sum(Shifted_profile.^2,1);
mean_sum_profile=mean(sum_profile);

variance_no_noise=zeros(1,number_of_columns);
for q = 1:1: number_of_columns
    sum_DS=sum((Shifted_profile(:,q).^2),1);
    ampPower(q)=sum_DS;
    if (sum_DS> mean_sum_profile)% first condition to satisfy dominant scatterer criteria
        Variance_remove_noise= var_rangebins(q);
        variance_no_noise(q)=Variance_remove_noise;
   
    end 
end
%figure;
%plot(ampPower);
%hold on
%yline(mean_sum_profile);
%hold off
%xlabel('Range(m)','fontsize',15)
%ylabel('Amplitude power','fontsize',15)
%set(gcf,'color','w')

variance_no_noise_no_zeros=nonzeros(variance_no_noise);

if isempty(variance_no_noise_no_zeros) %looks to see if range bins satisfy dominant scatterer criteria
    IsarImage=zeros(number_of_rows,number_of_columns);
    figure;
    Plot_IsarImage=imagesc(xAxis,yAxis,IsarImage);
    xlabel('Doppler frequency (Hz)')
    ylabel('Range(m)')
    axis xy;
    colormap('jet');
    colorbar
    axis ij;
    set(gcf,'color','w')

else 
    min_Variance=min(variance_no_noise_no_zeros); %second condition to satisfy dominant scatterer criteria
    pos_min_variance=find(var_rangebins==min_Variance);


    phase_aligned= angle(Aligned_HRR_profiles(:,pos_min_variance));
    %figure;
    %plot(phase_aligned); %plot phase history of range profile satisfying dominant scatter criteria
    %set(gcf,'color','w')
    %xlabel('Range(m)','fontsize',15)
    %ylabel('Phase history','fontsize',15)
    
    DS_CC= exp(-1i*phase_aligned);
    phaseComp_Matrix= repmat(DS_CC,1,number_of_columns); %deriving phase compensation matrix

    phase_Shift=Aligned_HRR_profiles.*phaseComp_Matrix;
    window_vector= hamming(number_of_rows);
    w=repmat(window_vector,1,number_of_columns);
    IsarImage=fftshift(fft(phase_Shift.*w,[],1),1);

    figure;
    Plot_IsarImage=imagesc(xAxis,yAxis,20*log10(abs(IsarImage))); %plot final ISAR image
    ylabel('Doppler frequency (Hz)','fontsize',15)
    xlabel('Range(m)','fontsize',15)
    axis xy;
    colormap('jet');
    colorbar
    set(gcf,'color','w')
    
end 
end 

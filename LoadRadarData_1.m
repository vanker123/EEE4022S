
clear all;
close all;

load('StrikemasterDataFormatted.mat');

HRRProfiles = Settings.HRR.HRR_calib_velcomppeak_win.'; 
[NumProfiles NumRangeBins] = size(HRRProfiles);

EffectivePRF = 1/Settings.IQData.timePatterns(2);
Bandwidth = (Settings.Pattern.FStep(1))*1e6*NumRangeBins;
RangeResolution = 3e8/(2*Bandwidth);

% Plot the HRR profiles 
RangeAxis_m = (0:1:(NumRangeBins-1))*RangeResolution;
ProfilesAxis = 1:1:NumProfiles;
figure; imagesc(RangeAxis_m,ProfilesAxis, 20*log10(abs(HRRProfiles)));
set(gca,'FontSize',14)
set(gcf,'color','w')
ylabel('Profile Number','fontsize',14);
xlabel('Range (m)','fontsize',14);
title('HRR profiles','fontsize',14);
colormap('jet');
colorbar;

%% Using window function to split dataset into smaller batches of results 

WindowLength = 64;
OverlapFactor=0.5;
overlap=0.5*WindowLength;

%subsetProfiles= HRRProfiles(1:561,:);      %results not too clear
%subsetProfiles= HRRProfiles(561:1122,:);   %results not good 
%subsetProfiles= HRRProfiles(1122:1683,:);  %results not good 
subsetProfiles= HRRProfiles(1683:2244,:);   %a few good results
%subsetProfiles= HRRProfiles(2244:2805,:);  %very bad results

y_length=size(subsetProfiles,1);
ShiftNextFrame = round(WindowLength*OverlapFactor);      
NumberofFrames = round((y_length - WindowLength)/ShiftNextFrame + 1); % Formula to calculate the total number of frames 

contrast= zeros(1,NumberofFrames-1);
entropy = zeros(1,NumberofFrames-1);
%%
for f = 1:NumberofFrames-1
    StartIdx = 1 + ShiftNextFrame*(f-1);
    StopIdx = 1+(WindowLength-1) + ShiftNextFrame*(f-1);
    subset= subsetProfiles(StartIdx:StopIdx,:);
    aligned_profiles= aligned_range(subset);

    [IsarImage Plot_ISAR]=calculate_autofocus(aligned_profiles,1:size(subset,1),(-WindowLength/2:1:(WindowLength/2-1))*EffectivePRF/WindowLength);
    
    Contrast_Value=calculate_contrast(Plot_ISAR);
    Entropy_Value=Entropy_of_ISARimage(Plot_ISAR);
    contrast(f)=Contrast_Value;
    entropy(f)=Entropy_Value;
    
    
end 
%maxc=find(contrast==max(contrast));
%mine=find(entropy==min(entropy));

close all
clear all
load('C:\Users\Abhishek Mukherjee\Downloads\DataIndex.mat');
load('C:\Users\Abhishek Mukherjee\Downloads\CData.mat');
load('C:\Users\Abhishek Mukherjee\Downloads\DataIndex5.mat');
load('C:\Users\Abhishek Mukherjee\Downloads\DataIndex7.mat');
load cdata.mat
load dataindex.mat
load DataIndex5.mat
load DataIndex7.mat

OriginalData = {user1; user2; user3; user5; user6; user7; user8; user9};
PeakIndices = {DataIndex1; DataIndex2; DataIndex3; DataIndex5; DataIndex6; DataIndex7; DataIndex8; DataIndex9};

Fs_original = 32;
t_upsample = 0:1/Fs_original:(length(user1)-1)/Fs_original;
Fs_downsample = 0.2;
% t_downsample = 0:1/Fs_downsample:((length(user1)/(Fs_original/Fs_downsample))-1)/Fs_downsample;

NumUsers = size(OriginalData,1);

Reconstruction = cell(NumUsers,1);
Splines = cell(NumUsers,1);
SplinePeaks = cell(NumUsers,1);
DownSample = cell(NumUsers,1);
FP = cell(NumUsers,1);

w = waitbar(0,'Initializing waitbar...');

for i=1:NumUsers
    perc = fix((i-1)/NumUsers*100);
    waitbar(perc/100,w,sprintf('Reconstructing %.0f of %.0f',i,NumUsers));
    Reconstruction{i} = ReconstructingSCR_PeakDetection(OriginalData{i},Fs_original,Fs_downsample,i);
    DownSample{i} = downsample(OriginalData{i},Fs_original/Fs_downsample);
    t_downsample = 0:1/Fs_downsample:(length(DownSample{i})-1)/Fs_downsample;
    Splines{i} = pchip(t_downsample,DownSample{i},t_upsample);
%     [~,SplinePeaks{i}] = findpeaks(Splines{i});
end

close(w)

[MSEReconstruction,NMAEReconstruction,NumMissedPeaksReconstuction] = PeakMSE(PeakIndices,Reconstruction,Fs_original,Fs_downsample);
[MSESplines,NMAESplines,NumMissedPeaksSplines] = PeakMSE(PeakIndices,Splines,Fs_original,Fs_downsample);
[MSEPD,NMAEPD,NumMissedPeaksPD] = PeakMSE(PeakIndices,DownSample,Fs_original,Fs_downsample,'DownSampled');

figure(1)
c = distinguishable_colors(3);
subplot(3,1,1)
b = bar([1;2;3;5;6;7;8;9],[[MSEReconstruction{:}]' [MSESplines{:}]' [MSEPD{:}]']);
for i=1:size(c,1)
    b(i).FaceColor = c(i,:);
end
title('MSE Between Labeled Peak Time and Calculated Peak Time')
legend('Proposed Method', 'PCHIP', 'PeakDetection')
xlabel('User')
ylabel('MSE (s^2)')
set(gca,'FontSize',16)
set(gca,'FontWeight','bold')
% set(gca,'YScale','log')
grid on
set(gca,'Gridalpha',0.6)

subplot(3,1,2)
b = bar([1;2;3;5;6;7;8;9],[[NMAEReconstruction{:}]' [NMAESplines{:}]' [NMAEPD{:}]']);
for i=1:size(c,1)
    b(i).FaceColor = c(i,:);
end
title('NMAE Between Labeled Peak Time and Calculated Peak Time')
legend('Proposed Method', 'PCHIP', 'PeakDetection')
xlabel('User')
ylabel('NMAE (samples at 0.2 Hz)')
set(gca,'FontSize',16)
set(gca,'FontWeight','bold')
% set(gca,'YScale','log')
grid on
set(gca,'Gridalpha',0.6)

subplot(3,1,3)
b = bar([1;2;3;5;6;7;8;9],[[NumMissedPeaksReconstuction{:}]' [NumMissedPeaksSplines{:}]' [NumMissedPeaksPD{:}]']);
for i=1:size(c,1)
    b(i).FaceColor = c(i,:);
end
title('Difference Between Number of Labeled Peaks and Calculated Peaks')
legend('Proposed Method', 'PCHIP', 'Peak Detection')
xlabel('User')
ylabel('Number of Missed Peaks')
set(gca,'FontSize',16)
set(gca,'FontWeight','bold')
% set(gca,'YScale','log')
grid on
set(gca,'Gridalpha',0.6)


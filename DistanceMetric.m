clear all
close all
load('C:\Users\Abhishek Mukherjee\Downloads\FoundPeaks');
load('C:\Users\Abhishek Mukherjee\Downloads\CData.mat');
load('C:\Users\Abhishek Mukherjee\Downloads\DataIndex.mat');
load('C:\Users\Abhishek Mukherjee\Downloads\CData.mat');
load('C:\Users\Abhishek Mukherjee\Downloads\DataIndex5.mat');
load('C:\Users\Abhishek Mukherjee\Downloads\DataIndex7.mat');
DataIndex={(DataIndex1),(DataIndex2),(DataIndex3),(DataIndex5),(DataIndex6),(DataIndex7),(DataIndex8),(DataIndex9)}
%% Proposed method:
SamplingFreqn=32;
for Num2=1:size((DataIndex),2)
DataIndexPrime=DataIndex{Num2};
PeakSampleInstants{Num2} = sort(DataIndexPrime(2:2:length(DataIndexPrime)));  
temp_RealLocations_peak = PeakSampleInstants{Num2};
DetectedLocations_peaks{Num2} = FoundPeaksReconstruction{Num2};
DetectedLocations_peaks{Num2} = EstimatedOnsetTimeInstants{Num2};
temp_DetectedLocations_peak = DetectedLocations_peaks{Num2};
for j=1:length(temp_DetectedLocations_peak)
 diff_in_time_peak_precise{1,j} = (abs(temp_DetectedLocations_peak(j) - temp_RealLocations_peak(:)))/SamplingFreqn;
 PrecisionMSE_proposed{Num2,j}=min(diff_in_time_peak_precise{1,j});
 end
for j=1:length(temp_RealLocations_peak)
 diff_in_time_peak_recall{1,j} = (abs(temp_RealLocations_peak(j) - temp_DetectedLocations_peak(:)))/SamplingFreqn;
 RecallMSE_proposed{Num2,j}=min(diff_in_time_peak_recall{1,j});
end
TotalPreci_proposed(Num2)=sum(cell2mat(PrecisionMSE_proposed(Num2,:)));
TotalRecall_proposed(Num2)=sum(cell2mat(RecallMSE_proposed(Num2,:)));
FscoresMSE_proposed(Num2)= (TotalPreci_proposed(Num2)*TotalRecall_proposed(Num2))/(TotalPreci_proposed(Num2)+TotalRecall_proposed(Num2));
end
%% splines
for Num2=1:size((DataIndex),2)
DataIndexPrime=DataIndex{Num2};
PeakSampleInstants{Num2} = sort(DataIndexPrime(1:2:length(DataIndexPrime)-1));  
temp_RealLocations_peak = PeakSampleInstants{Num2};
DetectedLocations_peaks{Num2} = FoundPeaksSplines{Num2};
temp_DetectedLocations_peak = DetectedLocations_peaks{Num2};
for j=1:length(temp_DetectedLocations_peak)
 diff_in_time_peak_precise{2,j} = (abs(temp_DetectedLocations_peak(j) - temp_RealLocations_peak(:)))/SamplingFreqn;
 PrecisionMSE_splines{Num2,j}=min(diff_in_time_peak_precise{2,j});
 end
for j=1:length(temp_RealLocations_peak)
 diff_in_time_peak_recall{2,j} = (abs(temp_RealLocations_peak(j) - temp_DetectedLocations_peak(:)))/SamplingFreqn;
 RecallMSE_splines{Num2,j}=min(diff_in_time_peak_recall{2,j});
end
TotalPreci_splines(Num2)=sum(cell2mat(PrecisionMSE_splines(Num2,:)));
TotalRecall_splines(Num2)=sum(cell2mat(RecallMSE_splines(Num2,:)));
FscoresMSE_splines(Num2)= (TotalPreci_splines(Num2)*TotalRecall_splines(Num2))/(TotalPreci_splines(Num2)+TotalRecall_splines(Num2));
end
%% PeakDetection
for Num2=1:size((DataIndex),2)
DataIndexPrime=DataIndex{Num2};
PeakSampleInstants{Num2} = sort(DataIndexPrime(1:2:length(DataIndexPrime)-1));  
temp_RealLocations_peak = PeakSampleInstants{Num2};
DetectedLocations_peaks{Num2} = FoundPeaksPD{Num2};
temp_DetectedLocations_peak = DetectedLocations_peaks{Num2};
for j=1:length(temp_DetectedLocations_peak)
 diff_in_time_peak_precise{3,j} = (abs(temp_DetectedLocations_peak(j) - temp_RealLocations_peak(:)))/SamplingFreqn;
 PrecisionMSE_PD{Num2,j}=min(diff_in_time_peak_precise{3,j});
 end
for j=1:length(temp_RealLocations_peak)
 diff_in_time_peak_recall{3,j} = (abs(temp_RealLocations_peak(j) - temp_DetectedLocations_peak(:)))/SamplingFreqn;
 RecallMSE_PD{Num2,j}=min(diff_in_time_peak_recall{3,j});
end
TotalPreci_PD(Num2)=sum(cell2mat(PrecisionMSE_PD(Num2,:)));
TotalRecall_PD(Num2)=sum(cell2mat(RecallMSE_PD(Num2,:)));
FscoresMSE_PD(Num2)= (TotalPreci_PD(Num2)*TotalRecall_PD(Num2))/(TotalPreci_PD(Num2)+TotalRecall_PD(Num2));
end
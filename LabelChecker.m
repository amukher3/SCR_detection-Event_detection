%% This code spits out a visual representation of the onset and the peaks 
% the data. One can use this to check the fidelity of the detected onsets
% and the peaks in the reconstruction. Just pass the detected onsets and
% the peaks in the DataIndex and the reconstruction in the user column. 
% This code assumes that the starting index is onset and the last label is
% a peak. 
%clear all
close all
% load('C:\Users\Abhishek Mukherjee\Downloads\DataIndex.mat');
% load('C:\Users\Abhishek Mukherjee\Downloads\CData.mat');
% load('C:\Users\Abhishek Mukherjee\Downloads\DataIndex7.mat');
% load('C:\Users\Abhishek Mukherjee\Downloads\DataIndex5.mat');
% load('C:\Users\Abhishek Mukherjee\Downloads\DataIndex4.mat');
% load('C:\Users\Abhishek Mukherjee\Downloads\user4_new.mat');
% load('C:\Users\Abhishek Mukherjee\Downloads\DataIndex4_incomplete');
%load('C:\Users\Abhishek Mukherjee\Downloads\DataIndex_Natrajan');
% load('C:\Users\Abhishek Mukherjee\Desktop\CompleteData');
% load('C:\Users\Abhishek Mukherjee\Desktop\DataIndex_s10');
% Data=CompleteData{9};

user= CroppedData;
DataIndex=MIT_DataIndex_167e4;

OnsetInstants=sort(DataIndex(1:2:length(DataIndex)-1));

for i=1:length(OnsetInstants)
epsilon(OnsetInstants(i))=user(OnsetInstants(i));
end

PeakInstants=sort(DataIndex(2:2:length(DataIndex)));

for i=1:length(PeakInstants)
oopsilon(PeakInstants(i))=user(PeakInstants(i));
end

figure;
plot(user)
hold on
plot(epsilon,'*')
figure
plot(user)
hold on
plot(oopsilon,'*')

%% Another method to check for false onsets and peaks in the reconstruction.

% One way to to locate a false onsets or a peak would be to rememeber the
% fact that there cannot be two consecutive onsets without a peak between
% them.Similar idea applies for peaks as well. 
% Hence, if we get two consecutive peaks or onsets we
% need to look for the maximam between them and just take that as our onset or peak.
% Onset is just the maximan of negative of the reconstruction. 
% If there is a FalsePeak or an FalseOnset "Idx" and "Idx_prime" would have
% two columns for that particular row. Indicating that index having the
% problem. 

% [ReconstructionFinal] = ReconstructingSCR_PeakDetection(user,SamplingFreqn,0.2,8);
% [MagOnset,OnsetSampleInstantstemp1] = findpeaks(-ReconstructionFinal);
%  [MagPeak,PeakSampleInstantstemp1] = findpeaks(ReconstructionFinal);
% 
% %FirstSample should be a peak.Just for reducing complexity of the different
% %cases. 
% for i=1:length(OnsetSampleInstantstemp1)
% if(PeakSampleInstantstemp1(1)>OnsetSampleInstantstemp1(i))
%  OnsetSampleInstantstemp(:)=OnsetSampleInstantstemp1(i+1:length(OnsetSampleInstantstemp1));
%  break;
% else
%  OnsetSampleInstantstemp(:)=OnsetSampleInstantstemp1(:); 
%  break;
% end
% end
% %LastSample should be an Onset.
% i=length(PeakSampleInstantstemp1);
% while(i>=1)
% if(PeakSampleInstantstemp1(i)>OnsetSampleInstantstemp(length(OnsetSampleInstantstemp)))
%  PeakSampleInstantstemp(:)=PeakSampleInstantstemp1(1:i-1);
%  break;
% else
%  PeakSampleInstantstemp(:)=PeakSampleInstantstemp1(:); 
%  break;
% end
% i=i-1;
% end
% %% For False Onsets
%  for i=1:length(PeakSampleInstantstemp)-1
% Indx(i,:)=find(PeakSampleInstantstemp(i)<OnsetSampleInstantstemp(:) & OnsetSampleInstantstemp(:)<PeakSampleInstantstemp(i+1));
%  end
% %% For FalsePeaks
% for i=1:length(OnsetSampleInstantstemp)-1
% Indx_prime(i,:)=find(OnsetSampleInstantstemp(i)<PeakSampleInstantstemp(:) & PeakSampleInstantstemp(:)<OnsetSampleInstantstemp(i+1));
% end
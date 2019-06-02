clear all
close all
load('C:\Users\amukher3\AppData\Local\Temp\temp3019822108431482480tmp\DataIndex.mat');
load('C:\Users\amukher3\AppData\Local\Temp\temp6491059322419438116tmp\CData.mat');
load('C:\Users\amukher3\Downloads\DataIndex7');
load('C:\Users\amukher3\Downloads\DataIndex5');
load('C:\Users\amukher3\Downloads\DataIndex4');
Users= {(user1),(user2),(user3),(user4),(user5),(user6),(user7),(user8),(user9)};
DataIndexes= {(DataIndex1),(DataIndex2),(DataIndex3),(DataIndex4),(DataIndex5),(DataIndex6),(DataIndex7),(DataIndex8),(DataIndex9)};
for Num1=1:length(Users)
user=Users{(Num1)};
DataIndex=DataIndexes{Num1};
SamplingFreqn=32;
OnsetSampleInstants = sort(DataIndex(2:2:length(DataIndex))); % ground truth
PeakSampleInstants = sort(DataIndex(1:2:length(DataIndex)-1)); % ground truth
WindowDuration=15; %time in seconds representing the duration during which we expect to observe the event given the stimulus. 
eventUseInd = [478,508,523,568,610,637];
%% For original signal
OnsetTimeInstants = OnsetSampleInstants/SamplingFreqn;
PeakTimeInstants = PeakSampleInstants/SamplingFreqn;
for i=1:length(eventUseInd)
    diff=OnsetTimeInstants - eventUseInd(i);
    Idx_original(i)=min(find(diff>0)); %looking for the closest onset after the stimuli 
    if((OnsetTimeInstants(Idx_original(i))-eventUseInd(i))< WindowDuration)
        NearestOnset(i)=OnsetTimeInstants(Idx_original(i)); % the nearest onset after the stimulus within WindowDuration
    else
        NearestOnset(i)=NaN; % if there are no events within the windowduration
    end
end
NearestOnset(find(isnan(NearestOnset)))=[]; %getting rid of the NaN instances

for i=1:length(NearestOnset)
    diff = PeakTimeInstants - NearestOnset(i);
    Idx_original_prime(i) = min(find(diff>0));
    NearestPeak(i)=PeakTimeInstants(Idx_original_prime(i));    % looking for the peak corresponding to a particular onset   
end

% After we have got the peak for the onset corresponding to a given
% stimuli we look for the next nearest onset after that peak
%  which is named as PseudoOnset.
% It is seen as Pseudo since it doesn't belong  to any stimuli instances. 

for i=1:length(NearestPeak)
    diff = OnsetTimeInstants - NearestPeak(i);
    Idx_original_tilda(i) = min(find(diff>0)); %closest pseudo onset after the peak
    NearestPseudoOnset(i)= OnsetTimeInstants(Idx_original_tilda(i));       
end
 FallTimeOriginal{Num1} = NearestPseudoOnset - NearestPeak;
 RiseTimeOriginal{Num1} = NearestPeak- NearestOnset;

%% For Reconstruction

[ReconstructionFinal]=ReconstructingSCRUsingPeakDetection(user)
[MagOnset,OnsetSampleInstantstemp1] = findpeaks(-ReconstructionFinal,'MinPeakProminence',0.1); %looking for the onsets in the reconstruction.
OnsetReconTimeInstants = OnsetSampleInstantstemp1/SamplingFreqn; %in seconds
[MagPeak,PeakSampleInstantstemp1] = findpeaks(ReconstructionFinal,'MinPeakProminence',0.1); %looking for peaks in the reconstruction. 
PeakReconTimeInstants = PeakSampleInstantstemp1/SamplingFreqn; %in seconds. 
for i=1:length(eventUseInd)
    diff = OnsetReconTimeInstants - eventUseInd(i);
    Idx(i)=min(find(diff>0));
    if((OnsetReconTimeInstants(Idx(i))-eventUseInd(i))< WindowDuration) %looking for the closest onset after the stimuli.
        NearestOnset_prime(i)=OnsetReconTimeInstants(Idx(i));
        else
        NearestOnset_prime(i)=NaN;
    end
end
NearestOnset_prime(find(isnan(NearestOnset_prime)))=[];

for i=1:length(NearestOnset_prime)
    diff = PeakReconTimeInstants - NearestOnset_prime(i);
    Idx_prime(i) = min(find(diff>0));
    NearestPeak_prime(i)=PeakReconTimeInstants(Idx_prime(i));      
end

for i=1:length(NearestPeak_prime)
    diff = OnsetReconTimeInstants - NearestPeak_prime(i);
    Idx_tilda(i) = min(find(diff>0));
    NearestPseudoOnset_prime(i)= OnsetReconTimeInstants(Idx_tilda(i));      
end
 FallTimeRecon{Num1} = NearestPseudoOnset_prime - NearestPeak_prime;
 RiseTimeRecon{Num1} = NearestPeak_prime - NearestOnset_prime;
 [Pvalue_FallTime(Num1),h_FallTime(Num1)]=ranksum(FallTimeOriginal{Num1},FallTimeRecon{Num1});
 [Pvalue_RiseTime(Num1),h_RiseTime(Num1)]=ranksum(RiseTimeOriginal{Num1},RiseTimeRecon{Num1});
clearvars -except  FallTimeOriginal RiseTimeOriginal Users DataIndexes FallTimeRecon RiseTimeRecon Pvalue_FallTime h_FallTime Pvalue_RiseTime h_RiseTime

end
 y=[mean(RiseTimeOriginal{1,1}) mean(RiseTimeRecon{1,1});mean(RiseTimeOriginal{1,2}) mean(RiseTimeRecon{1,2});mean(RiseTimeOriginal{1,3}) mean(RiseTimeRecon{1,3});mean(RiseTimeOriginal{1,4}) mean(RiseTimeRecon{1,4});mean(RiseTimeOriginal{1,5}) mean(RiseTimeRecon{1,5});mean(RiseTimeOriginal{1,6}) mean(RiseTimeRecon{1,6});mean(RiseTimeOriginal{1,7}) mean(RiseTimeRecon{1,7});mean(RiseTimeOriginal{1,8}) mean(RiseTimeRecon{1,8});mean(RiseTimeOriginal{1,9}) mean(RiseTimeRecon{1,9});];
 figure;
 bar(y);
 xlabel('users');
 ylabel('Average Rise Time');
 legend('Original','Reconstruction');
 title('RiseTime')
 yBar=[mean(FallTimeOriginal{1,1}) mean(FallTimeRecon{1,1});mean(FallTimeOriginal{1,2}) mean(FallTimeRecon{1,2});mean(FallTimeOriginal{1,3}) mean(FallTimeRecon{1,3});mean(FallTimeOriginal{1,4}) mean(FallTimeRecon{1,4});mean(FallTimeOriginal{1,5}) mean(FallTimeRecon{1,5});mean(FallTimeOriginal{1,6}) mean(FallTimeRecon{1,6});mean(FallTimeOriginal{1,7}) mean(FallTimeRecon{1,7});mean(FallTimeOriginal{1,8}) mean(FallTimeRecon{1,8});mean(FallTimeOriginal{1,9}) mean(FallTimeRecon{1,9});];
 figure;
 bar(yBar);
 xlabel('users');
 ylabel('Average Fall Time');
 legend('Original','Reconstruction');
 title('FallTime')
%% One way ANOVA for RiseTime and FallTime

%There are 2 separate design matrices. FallTime and RiseTime.
%The design matrix has 2 columns signifying the 2 groups - Original and
%Reconstruction.

% ReconFallTimeVect=[mean(FallTimeRecon{1,1}),mean(FallTimeRecon{1,2}),mean(FallTimeRecon{1,3}),mean(FallTimeRecon{1,4}),mean(FallTimeRecon{1,5}),mean(FallTimeRecon{1,6}),mean(FallTimeRecon{1,7}),mean(FallTimeRecon{1,8})];
% OriginalFallTimeVect=[mean(FallTimeOriginal{1,1}),mean(FallTimeOriginal{1,2}),mean(FallTimeOriginal{1,3}),mean(FallTimeOriginal{1,4}),mean(FallTimeOriginal{1,5}),mean(FallTimeOriginal{1,6}),mean(FallTimeOriginal{1,7}),mean(FallTimeOriginal{1,8})];
% ReconFallTimeVect=[ReconFallTimeVect NaN(1,length(OriginalFallTimeVect)-length(ReconFallTimeVect))];
% DesignMat_FallTime=[OriginalFallTimeVect',ReconFallTimeVect'];
% [pValueFallTime]= anova1(DesignMat_FallTime)
% ReconRiseTimeVect=[mean(RiseTimeRecon{1,1}),mean(RiseTimeRecon{1,2}),mean(RiseTimeRecon{1,3}),mean(RiseTimeRecon{1,4}),mean(RiseTimeRecon{1,5}),mean(RiseTimeRecon{1,6}),mean(RiseTimeRecon{1,7}),mean(RiseTimeRecon{1,8})];
% OriginalRiseTimeVect=[mean(RiseTimeOriginal{1,1}),mean(RiseTimeOriginal{1,2}),mean(RiseTimeOriginal{1,3}),mean(RiseTimeOriginal{1,4}),mean(RiseTimeOriginal{1,5}),mean(RiseTimeOriginal{1,6}),mean(RiseTimeOriginal{1,7}),mean(RiseTimeOriginal{1,8})];
% ReconRiseTimeVect=[ReconRiseTimeVect NaN(1,length(OriginalRiseTimeVect)-length(ReconRiseTimeVect))];
% DesignMat_RiseTime=[OriginalRiseTimeVect',ReconRiseTimeVect'];
% [pValueRiseTime]= anova1(DesignMat_RiseTime)

%% Wilcoxon rank sum or Unbalanced ANOVA
% ReconFallTimeVect=[(FallTimeRecon{1,1}),(FallTimeRecon{1,2}),(FallTimeRecon{1,3}),(FallTimeRecon{1,4}),(FallTimeRecon{1,5}),(FallTimeRecon{1,6}),(FallTimeRecon{1,7}),(FallTimeRecon{1,8})];
% OriginalFallTimeVect=[(FallTimeOriginal{1,1}),(FallTimeOriginal{1,2}),(FallTimeOriginal{1,3}),(FallTimeOriginal{1,4}),(FallTimeOriginal{1,5}),(FallTimeOriginal{1,6}),(FallTimeOriginal{1,7}),(FallTimeOriginal{1,8})];
% ReconFallTimeVect=[ReconFallTimeVect NaN(1,length(OriginalFallTimeVect)-length(ReconFallTimeVect))];
% DesignMat_FallTime = [OriginalFallTimeVect',ReconFallTimeVect'];
% [pValueFallTime]= anova1(DesignMat_FallTime)
% ReconRiseTimeVect=[(RiseTimeRecon{1,1}),(RiseTimeRecon{1,2}),(RiseTimeRecon{1,3}),(RiseTimeRecon{1,4}),(RiseTimeRecon{1,5}),(RiseTimeRecon{1,6}),(RiseTimeRecon{1,7}),(RiseTimeRecon{1,8})];
% OriginalRiseTimeVect = [(RiseTimeOriginal{1,1}),(RiseTimeOriginal{1,2}),(RiseTimeOriginal{1,3}),(RiseTimeOriginal{1,4}),(RiseTimeOriginal{1,5}),(RiseTimeOriginal{1,6}),(RiseTimeOriginal{1,7}),(RiseTimeOriginal{1,8})];
% ReconRiseTimeVect = [ReconRiseTimeVect NaN(1,length(OriginalRiseTimeVect)-length(ReconRiseTimeVect))];
% DesignMat_RiseTime = [OriginalRiseTimeVect',ReconRiseTimeVect'];
% [pValueRiseTime]= anova1(DesignMat_RiseTime)



function [MSE,NormalizedMeanAbsoluteError,NumMissedPeaks] = PeakMSE(PeakIndices,Data,Fs_original,Fs_downsample,varargin)
% Calculates the MSE between the peak times of a reconstruction and
% expertly labeled peak indices. Also gives the number of missed peaks
% labeled in the orignal data but not found in the reconstruction.

FoundPeaks = cell(1,size(Data,1));
MSE = cell(1,size(Data,1));
NumMissedPeaks = cell(1,size(Data,1));
NormalizedMeanAbsoluteError = cell(1,size(Data,1));

for i=1:size(PeakIndices,1)
    PeakIndices{i} = PeakIndices{i}(1:2:size(PeakIndices{i},1)-1);
    [~,FoundPeaks{i}] = findpeaks(Data{i});
    if ~isempty(varargin)
        if varargin{1} == 'DownSampled'
            FoundPeaks{i} = FoundPeaks{i} *(Fs_original/Fs_downsample);
        end
    end 
    MSE{i} = 0;
    NormalizedMeanAbsoluteError{i} = 0;
    if size(FoundPeaks{i},1) < size(PeakIndices{i},1)
        for j=1:size(FoundPeaks{i},1)
            [MinDistance] = min(abs(PeakIndices{i}-FoundPeaks{i}(j)));
            MSE{i} = MSE{i}+(MinDistance/Fs_original)^2;
            NormalizedMeanAbsoluteError{i} = NormalizedMeanAbsoluteError{i} + abs(MinDistance/Fs_original);
        end
        MSE{i} = MSE{i}/(size(FoundPeaks{i},1));
        NormalizedMeanAbsoluteError{i} = Fs_downsample*NormalizedMeanAbsoluteError{i}/(size(FoundPeaks{i},1));
    else
        for j=1:size(PeakIndices{i},1)
            [MinDistance] = min(abs(FoundPeaks{i}-PeakIndices{i}(j)));
            MSE{i} = MSE{i}+(MinDistance/Fs_original)^2;
            NormalizedMeanAbsoluteError{i} = NormalizedMeanAbsoluteError{i} + abs(MinDistance/Fs_original);
        end
        MSE{i} = MSE{i}/(size(PeakIndices{i},1));
        NormalizedMeanAbsoluteError{i} = Fs_downsample*NormalizedMeanAbsoluteError{i}/(size(FoundPeaks{i},1));
    end
    
    NumMissedPeaks{i} = size(PeakIndices{i},1) - size(FoundPeaks{i},1); 
end


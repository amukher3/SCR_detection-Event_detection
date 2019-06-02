function halfRecTimes = computeHalfRecTimes(sc,freq,peakLocs,onsetLocs, ...
    amplitudes)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

nPeaks = length(peakLocs);
% Compute half-recovery times of SCRs where possible. If the first instance
% of a half-recovery SC measurement is after the start of another SCR, then
% treat the half-recovery time as NaN
halfRecSc = sc(peakLocs) - amplitudes/2;
possibleHalfRecLocs = zeros(1,nPeaks);
halfRecTimes = zeros(1,nPeaks);
for iPeak = 1:nPeaks
    % Find first sample after specified peak that drops below the
    % half-recovery SC measurement. If there is none (e.g. end of data set),
    % then set half-recovery time to NaN.
    halfRecLocAfterPeak = find(sc(peakLocs(iPeak):end) < halfRecSc(iPeak),1);
    if isempty(halfRecLocAfterPeak)
        halfRecTimes(iPeak) = NaN;
        continue
    end    
    possibleHalfRecLocs(iPeak) = halfRecLocAfterPeak + peakLocs(iPeak) - 1;
    if iPeak < nPeaks
        % Check if possible half-recovery time occurs beyond next SCR;
        % if so, then discard it because of overlapping SCR
        if possibleHalfRecLocs(iPeak) > onsetLocs(iPeak+1)
            halfRecTimes(iPeak) = NaN;
            continue
        end
    end
    
    halfRecTimes(iPeak) = (possibleHalfRecLocs(iPeak)-peakLocs(iPeak)) / freq;
end

end


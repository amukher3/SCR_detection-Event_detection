%function [obs,est,matchLocs] = computeScrFeatures(scObs,scEst,t, ...
%    obsOnsetTimes,obsPeakTimes)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% clear all
% close all
% load('C:\Users\Abhishek Mukherjee\Downloads\DataIndex.mat');
% load('C:\Users\Abhishek Mukherjee\Downloads\CData.mat');
% load('C:\Users\Abhishek Mukherjee\Downloads\DataIndex5.mat');
% load('C:\Users\Abhishek Mukherjee\Downloads\DataIndex4.mat');
% load('C:\Users\Abhishek Mukherjee\Downloads\DataIndex7.mat');
function [obs,est,matchLocs] = computeScrFeatures(scObs,DataIndex)
% scObs=user2;
% DataIndex=DataIndex2;
freq = 32;
[scEst]= ReconstructingSCR_PeakDetection(scObs)
[~,estOnset] = findpeaks(-scEst,'MinPeakProminence',.01);
estOnsetTimes = estOnset/freq;
[~,estPeak] = findpeaks(scEst,'MinPeakProminence',.01);
estPeakTimes = estPeak/freq;
minLength = min([length(estOnsetTimes)],[length(estPeakTimes)]);
estOnsetTimes=estOnsetTimes(1:minLength);
estPeakTimes=estPeakTimes(1:minLength);
obsOnsetTimes=((sort(DataIndex(2:2:length(DataIndex))))/freq);
obsPeakTimes=((sort(DataIndex(1:2:length(DataIndex)-1)))/freq);
%t = 0:1/freq:duration-1/freq;
%if nargin < 8
    maxDist = 3;
%end
obs.peakTimes = obsPeakTimes;
obs.peakLocs = round(obs.peakTimes*freq)+1;
obs.onsetTimes = obsOnsetTimes;
obs.onsetLocs = round(obs.onsetTimes*freq)+1;
est.onsetTimes = estOnsetTimes;
est.onsetLocs = round(est.onsetTimes*freq)+1;
est.peakTimes = estPeakTimes;
est.peakLocs = round(est.peakTimes*freq)+1;

% Filter observed (noisy) skin conductance using median filter
% scObsFilt = medfilt1(scObs,3);
scObsFilt = scObs;

% Match up each estimated peak to the closest observed peak
[estPeakTimesGrid,obsPeakTimesGrid] = ndgrid(est.peakTimes,obs.peakTimes);
estPeakDist = abs(estPeakTimesGrid - obsPeakTimesGrid);
[estPeakMatchDist,matchLocs] = min(estPeakDist,[],2);
peaksOutsideMaxDist = find(estPeakMatchDist > maxDist);
if ~isempty(peaksOutsideMaxDist)
    fprintf(['Estimated peaks at the following times have match distance ' ...
        '> %4.2f seconds:\n'],maxDist);
    fprintf('%4.2f\n',est.peakTimes(peaksOutsideMaxDist));
end

% Compute SCR rise times
obs.riseTimes = obs.peakTimes - obs.onsetTimes;
est.riseTimes = est.peakTimes - est.onsetTimes;

% Compute SCR amplitudes
obs.amplitudes = scObsFilt(obs.peakLocs) - scObsFilt(obs.onsetLocs);
est.amplitudes = scEst(est.peakLocs) - scEst(est.onsetLocs);

% Compute half-recovery times of SCRs where possible. If the first instance
% of a half-recovery SC measurement is after the start of another SCR, then
% treat the half-recovery time as NaN
obs.halfRecTimes = computeHalfRecTimes(scObsFilt,freq,obs.peakLocs, ...
    obs.onsetLocs,obs.amplitudes);
est.halfRecTimes = computeHalfRecTimes(scEst,freq,est.peakLocs, ...
    est.onsetLocs,est.amplitudes);
end





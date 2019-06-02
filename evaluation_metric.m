%This function calculates the Precision and Recall with respect to a
%distance parameter from the event instats i.e Precision/Recall vs.
%distance. 
clear all
close all
load('C:\Users\Abhishek Mukherjee\Downloads\DataIndex.mat');
load('C:\Users\Abhishek Mukherjee\Downloads\CData.mat');
load('C:\Users\Abhishek Mukherjee\Downloads\DataIndex5.mat');
% Users= {(user1),(user2),(user3),(user5),(user6),(user8),(user9)};
%DataIndexes= {(DataIndex1),(DataIndex2),(DataIndex3),(DataIndex5),(DataIndex6),(DataIndex8),(DataIndex9)};
 Users = {(user2)};
 DataIndexes = {(DataIndex2)};
%% Evaluation metrics for the proposed algorithm. 
S_d=0;S_r=0;
 for Num2=1:length(Users)
%[Likelihood_under_Null,AUC,X,Y,scores,Threshold,Likelihood_under_Alternative,mle_estimates,Window,NoiseWindows,EventWindows,condn,OnsetTimeInstants,tau_mle,SamplingFreqn,Delta,T,OnsetSampleInstants,WindowTimeLength,labels,noise_var,TimeInstantOfFirstSample,locs] = DevPeakLocations(Users{:,Num2},DataIndexes{:,Num2});
[Likelihood_under_Null,AUC,scores,Likelihood_under_Alternative,Window,OnsetTimeInstants,tau_mle,SamplingFreqn,Delta,T,OnsetSampleInstants,WindowTimeLength,labels,noise_var,TimeInstantOfFirstSample,locs,mle_estimates] = DevPeakLocations(Users{:,Num2},DataIndexes{:,Num2});
estimates{Num2} = mle_estimates(:,:); % MLE estimates for the TimeWindows for every user. 
User_AUC(1,Num2)=AUC; % AUC's for the different users.
UserScores(:,Num2) = scores;
UserWindows{Num2}=Window;
Likelihood_under_Alternative_users(:,Num2) = Likelihood_under_Alternative(:);
Likelihood_under_Null_users(:,Num2) = Likelihood_under_Null(:);
DownsamplingFactor = T*SamplingFreqn;
locations{Num2} = locs; %locations of the peaks in the dev. stats.
Tau_MLE_locs{Num2} = tau_mle(locations{Num2}); % values of tau for the windows where the dev. statistics peaks. 
%DetectedLocations{Num2} = (Delta+(locations{Num2}*DownsamplingFactor-1)*(1/SamplingFreqn))- (Tau_MLE_locs{Num2}); %calculated at high sample rate
DetectedLocations{Num2} = (Delta+(locations{Num2})*T)-(Tau_MLE_locs{Num2}); % for low sampling rate
RealLocations{Num2} = OnsetTimeInstants(:);
RealLocationsSamples{Num2} = OnsetSampleInstants(:); % Onset Sample Instants(Handlabelled)
S_d=S_d+length(DetectedLocations{Num2});
S_r=S_r+length(RealLocations{Num2});  
end
dmax=1:1:5; % vary this to see performance at larger distances. 
Sdr=zeros(1,length(dmax)); 
Srd=zeros(1,length(dmax)); 
for i=1:length(dmax)
for Num2=1:length(Users)
temp_RealLocations = RealLocations{Num2};
temp_DetectedLocations = DetectedLocations{Num2};
for j=1:length(temp_DetectedLocations)
 diff_in_time{j} = temp_DetectedLocations(j) - temp_RealLocations(:);
 if (any(abs(diff_in_time{j})<=dmax(i)))
 Sdr(i) = Sdr(i)+1;
 end
end
for j=1:length(temp_RealLocations)
 diff_in_time{j} = temp_RealLocations(j) - temp_DetectedLocations(:);
 if (any(abs(diff_in_time{j})<=dmax(i)))
 Srd(i)=Srd(i)+1;
 end
end
end
end
preci=Sdr/S_d;
recall=Srd/S_r;
Fscore(:)=2*(preci(:).*recall(:))./(preci(:)+recall(:));

%% Evaluation metric using the peaks:
S_d_peak=0; S_r_peak=0;
Srd_peak=zeros(1,length(dmax)); 
Sdr_peak=zeros(1,length(dmax)); 
DownsamplingFactor = T*SamplingFreqn;
for Num2=1:length(Users)
user(:)= Users{Num2};
if(Delta==0)
DownUser = downsample(user,DownsamplingFactor);
else
DownUser = user(round(Delta*SamplingFreqn):DownsamplingFactor:length(user));
end
DownSampledUsers{Num2} = DownUser;
%UpSampledUser{Num2,:}=interp(DownUser,DownsamplingFactor); %For interpolated downsampled data.
[~,locPeakDetection]= findpeaks(DownSampledUsers{Num2});
DetectedPeaks{Num2}= locPeakDetection;
%[~,locPeakDetection] = findpeaks(UpSampledUser{Num2,:}); 
DetectedLocations_peak{Num2} = ((locPeakDetection-1)/(1/T));
DataIndex = DataIndexes{Num2};
PeakSamplingInstants{Num2} = sort(DataIndex(1:2:length(DataIndex)-1)); % comparing with the labeled peaks
PeakSamplingInstantsTemps = PeakSamplingInstants{Num2}; 
PeaksTimeInstants{Num2} =((PeakSamplingInstantsTemps-1)/SamplingFreqn);
RealLocations_peak{Num2} = PeaksTimeInstants{Num2};
S_d_peak=S_d_peak+length(DetectedLocations_peak{Num2});
S_r_peak= S_r_peak+length(RealLocations_peak{Num2});
end
for i=1:length(dmax)
for Num2=1:length(Users)
temp_RealLocations_peak = RealLocations_peak{Num2};
temp_DetectedLocations_peak = DetectedLocations_peak{Num2};
for j=1:length(temp_DetectedLocations_peak)
 diff_in_time_peak{j} = temp_DetectedLocations_peak(j) - temp_RealLocations_peak(:);
 if (any(abs(diff_in_time_peak{j})<=dmax(i))) 
 Sdr_peak(i)=Sdr_peak(i)+1;
 end
end
for j=1:length(temp_RealLocations_peak)
 diff_in_time_peak{j} = temp_RealLocations_peak(j) - temp_DetectedLocations_peak(:);
 if (any(abs(diff_in_time_peak{j})<=dmax(i)))
 Srd_peak(i)=Srd_peak(i)+1;
 end
end
end
end
preci_peak_detection=Sdr_peak/S_d_peak;
recall_peak_detection=Srd_peak/S_r_peak;
Fscore_peak_detection(:)=2*(preci_peak_detection(:).*recall_peak_detection(:))./(preci_peak_detection(:)+recall_peak_detection(:));
figure;
plot(dmax,preci,'--o');
hold on
plot(dmax,preci_peak_detection,'--o');
xlabel('dmax in sec');
ylabel('Precision');
legend('Proposed Method','Peak Detection');
figure;
plot(dmax,recall,'--o');
hold on
plot(dmax,recall_peak_detection,'--o');
xlabel('dmax in sec');
ylabel('Recall');
legend('Proposed Method','Peak Detection');
figure;
plot(dmax,Fscore,'--o');
hold on
plot(dmax,Fscore_peak_detection,'--o');
xlabel('dmax in sec');
ylabel('Fscore');
legend('Proposed Method','Peak Detection');
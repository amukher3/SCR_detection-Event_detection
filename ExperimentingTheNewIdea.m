%This function calculates the Precision and Recall with respect to a
%distance parameter from the event instats i.e Precision/Recall vs.
%distance. 
clear all
close all
load('C:\Users\Abhishek Mukherjee\Downloads\DataIndex.mat');
load('C:\Users\Abhishek Mukherjee\Downloads\CData.mat');
load('C:\Users\Abhishek Mukherjee\Downloads\DataIndex5.mat');
load('C:\Users\Abhishek Mukherjee\Downloads\DataIndex7.mat');
Users= {(user1),(user2),(user3),(user5),(user6),(user7),(user8),(user9)};
DataIndexes= {(DataIndex1),(DataIndex2),(DataIndex3),(DataIndex5),(DataIndex6),(DataIndex7),(DataIndex8),(DataIndex9)};
%Users= {(user1)};
%DataIndexes= {(DataIndex1)};
%% Detecting the peaks in the downsampled data:
T=4;
SamplingFreqn=32;
Delta=0; % vary delta to see the performance with shift in the experiments. 
DownsamplingFactor = T*SamplingFreqn;
for Num2=1:length(Users)
user= Users{Num2};
if(Delta==0)
DownUser = downsample(user,DownsamplingFactor);
else
DownUser = user(round(Delta*SamplingFreqn-1):DownsamplingFactor:length(user));
end
DownSampledUsers{Num2} = DownUser;
%UpSampledUser{Num2,:}=interp(DownUser,DownsamplingFactor); %For interpolated downsampled data.
[~,locPeakDetection]= findpeaks(DownSampledUsers{Num2});
%NumOfComponents=2;
% Y=fft(DownUser);
% SclFreqDomain = Y(1:NumOfComponents);
% SclTimeDomain=abs(ifft(SclFreqDomain,length(DownUser)));
% %SCLinTime=interp(SclTimeDomain,DownsamplingFactor);
% DownSampledUsers{Num2}=DownSampledUsers{Num2}-SclTimeDomain+10;
%DownSampledUsers{Num2}=DownSampledUsers{Num2}+100;
OnsetSample = locPeakDetection-1; % due to downsampling the onset sample is the detected peak in the downsampled data. 
FirstSample{Num2} = OnsetSample;
%[~,locPeakDetection] = findpeaks(UpSampledUser{Num2,:}); 
DataIndex = DataIndexes{Num2};
OnsetSampleInstants{Num2} = sort(DataIndex(2:2:length(DataIndex))); % even indexes are Onset time instants.The last Index is a Onset time instant and after that every alternate.
OnsetTimeInstants{Num2} = (OnsetSampleInstants{Num2}-1)/SamplingFreqn;
end

%% Making the windows
WindowSize=4; %Num of Samples in the window -1

 for i=1:length(Users)
     temp1 = FirstSample{i};
     DownUser = DownSampledUsers{i};
     WindowTemp=zeros(WindowSize+1,length(temp1));
 for j=1:length(temp1)
     if(temp1(j) + WindowSize > length(DownUser))
       PadSize = (temp1(j)+WindowSize) - length(DownUser);
       WindowTemp(:,j) = [DownUser(temp1(j):length(DownUser));zeros(PadSize,1)] ;
     else
       WindowTemp(:,j) = DownUser(temp1(j):temp1(j)+WindowSize) ;
     end
 end
  Window{i} = WindowTemp;
 end
 
%% Parameter estimation
A=[];
b=[];
lb= [.01,0.01]; % lower bound of the estimates:[alpha tau TauOne TauTwo]
ub= [24,T-.001]; %upper bound of the estimates
Aeq = []; %for equality constraint
beq = []; %for equality constraint
iter=10; % Num of independent runs
NumOfParameters=2;
TauOne=10; %  Time constants choosen empirically. 
TauTwo=1;
x=zeros(iter,NumOfParameters);
x0(1,:) = [lb(1)/2,.01]; % initial estiamtes of the parameters
x0(2,:) = [lb(1)/3,1];
x0(3,:) = [lb(1)/5,2];
x0(4,:) = [lb(1)/10,3];
x0(5,:) = [lb(1)/4,4];
x0(6,:) = [ub(1)/6,.01];
x0(7,:) = [ub(1)/ub(1),0.5];
x0(8,:) = [ub(1),1.5];
x0(9,:) = [ub(1),2.5];
x0(10,:) = [ub(1)/2,3.5];
x0(11,:) = [ub(1)/3,4.5];
x0(12,:) = [(lb(1)-ub(1))/2,1.2];
x0(13,:) = [(lb(1)-ub(1))/6,2.2];
x0(14,:) = [(lb(1)- ub(1))/8,3.2];
x0(15,:) = [(lb(1)-ub(1))/3,4.6];
S_d=0;S_r=0;
for Num=1:length(Users) %Num of users
    UserWindows=zeros((WindowSize+1),length(Window{1,Num}));
    UserWindows=Window{1,Num};
    NumOfWindows = size(UserWindows,2);
    alpha_mle=zeros(1,NumOfWindows);
    tau_mle=zeros(1,NumOfWindows);
for runs=1:NumOfWindows % Num of Windows in a user
  samples(:,runs) = UserWindows(:,runs);
    for i=1:iter 
options = optimoptions('fmincon','SpecifyObjectiveGradient',true); 
[x(i,:),val(runs,i),flag(runs,i)] = fmincon(@(t_prime)LeastSquareEstimates(t_prime,samples(:,runs),T,TauOne,TauTwo),x0(i,:),A,b,Aeq,beq,lb,ub,[],options);
    end
[~,I] = min(val(runs,:)); % checking for the minm. value of the Likelihood function
mle_estimates(runs,:)= x(I,:); % getting the estimate values(alpha,tau,tauOne,tauTwo) for the corresponding minm. value of the likelihood
alpha_mle(runs)= mle_estimates(runs,1);
tau_mle(runs) = mle_estimates(runs,2);
end
AlphaEstimatesUsers{Num} = alpha_mle;
TauEstimatesUsers{Num} = tau_mle;
EstimatedOnsetTimeInstants{Num} = ((FirstSample{Num})/(1/T)) - (TauEstimatesUsers{Num})';
S_d=S_d+length(EstimatedOnsetTimeInstants{Num});
S_r=S_r+length(OnsetTimeInstants{Num});  
end
  
%% Evaluation 
% Evaluation is done wrt the onset time instants for the proposed
% algorithm. Since we are trying to estimate delta/tau from the data.

dmax=1:1:4;  % vary this to see performance at larger distances. 
Sdr=zeros(1,length(dmax)); 
Srd=zeros(1,length(dmax)); 

for i=1:length(dmax)
for Num2=1:length(Users)
temp_RealLocations = OnsetTimeInstants{Num2};
temp_DetectedLocations = EstimatedOnsetTimeInstants{Num2};

for j=1:length(temp_DetectedLocations)
 diff_in_time{Num2,j} = temp_DetectedLocations(j) - temp_RealLocations(:);
 if (any(abs(diff_in_time{Num2,j})<=dmax(i)))
 Sdr(i) = Sdr(i)+1;
 end
end
for j=1:length(temp_RealLocations)
 diff_in_time{Num2,j} = temp_RealLocations(j) - temp_DetectedLocations(:);
 if (any(abs(diff_in_time{Num2,j})<=dmax(i)))
 Srd(i)=Srd(i)+1;
 end
end
end
end
preci=Sdr/S_d;
recall=Srd/S_r;
Fscore=2*(preci.*recall)./(preci+recall);
%% Evaluation metric using the peaks:
S_d_peak=0; S_r_peak=0;
Srd_peak=zeros(1,length(dmax)); 
Sdr_peak=zeros(1,length(dmax)); 
DownsamplingFactor = T*SamplingFreqn;
for Num2=1:length(Users)
user= Users{Num2};
%% comment this section while getting peaks from other methods.     
if(Delta==0)
DownUser = downsample(user,DownsamplingFactor);
else
DownUser = user(round(Delta*SamplingFreqn):DownsamplingFactor:length(user));
end
DownSampledUsers{Num2} = DownUser;
[~,locPeakDetection]= findpeaks(DownSampledUsers{Num2});
%% Use this part when getting the peak locations from other algorithms. 
DownSamplFreq = 1/T; % change this to get the locations at higher frequencies. 
DetectedPeaks{Num2}= locPeakDetection; % just throw in the detected peaks of other reconstructions here. 
DetectedLocations_peak{Num2} = ((locPeakDetection-1)/DownSamplFreq); % in seconds
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
 diff_in_time_peak{Num2,j} = temp_DetectedLocations_peak(j) - temp_RealLocations_peak(:);
 if (any(abs(diff_in_time_peak{Num2,j})<=dmax(i))) 
 Sdr_peak(i)=Sdr_peak(i)+1;
 end
end
for j=1:length(temp_RealLocations_peak)
 diff_in_time_peak{Num2,j} = temp_RealLocations_peak(j) - temp_DetectedLocations_peak(:);
 if (any(abs(diff_in_time_peak{Num2,j})<=dmax(i)))
 Srd_peak(i)=Srd_peak(i)+1;
 end
end
end
end
preci_peak_detection=Sdr_peak/S_d_peak;
recall_peak_detection=Srd_peak/S_r_peak;
Fscore_peak_detection=2*(preci_peak_detection.*recall_peak_detection)./(preci_peak_detection+recall_peak_detection);
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


% This function detections the peaks in the Likelihoods/difference of the
% Likelihoods(deviance statistics).
%Computes AUC by using the peaks in the deviance statistics as scores in
%perfcurve.
%function [Likelihood_under_Null,AUC,X,Y,scores,Threshold,Likelihood_under_Alternative,mle_estimates,Window,NoiseWindows,EventWindows,condn,OnsetTimeInstants,tau_mle,SamplingFreqn,Delta,T,OnsetSampleInstants,WindowTimeLength,labels,noise_var,TimeInstantOfFirstSample,loc] = DevPeakLocations(user,DataIndex)
function [Likelihood_under_Null,AUC,scores,Likelihood_under_Alternative,Window,OnsetTimeInstants,tau_mle,SamplingFreqn,Delta,T,OnsetSampleInstants,WindowTimeLength,labels,noise_var,TimeInstantOfFirstSample,loc,mle_estimates] = DevPeakLocations(user,DataIndex);
[Window,NoiseWindows,EventWindows,condn,OnsetTimeInstants,tau,SamplingFreqn,Delta,T,OnsetSampleInstants,WindowTimeLength,labels,TimeInstantOfFirstSample]= MakingTimeWindows(user,DataIndex);
TotalSizeOfWindow=size(Window,2); % Total number of windows 
NumOfParameters=4;
iter=10; % no. of independent initialisations
val=zeros(TotalSizeOfWindow,iter); %values of the Likelihood for the windows for each independent intialisation
flag=zeros(1,TotalSizeOfWindow);
NumOfSamplesInWindow = WindowTimeLength/T;
%sample_differences=zeros(1,NumOfSamplesInWindow-1); % when doing pairwise
%difference for the samples in the window. 
sample_differences=zeros(1,NumOfSamplesInWindow);
scores=zeros(1,TotalSizeOfWindow);
mle_estimates=zeros(TotalSizeOfWindow,NumOfParameters);
alpha_mle=zeros(1,TotalSizeOfWindow);
tau_mle=zeros(1,TotalSizeOfWindow);
tauOne_mle=zeros(1,TotalSizeOfWindow);
tauTwo_mle=zeros(1,TotalSizeOfWindow);
Likelihood_under_Null=zeros(1,TotalSizeOfWindow);
Likelihood_under_Alternative=zeros(1,TotalSizeOfWindow);
E=1;N=1; noise_var=1e-2;
for runs=1:TotalSizeOfWindow
%% Separating the windows into Noise windows and Event windows  
    if labels(runs)==0
        ys=NoiseWindows(:,N);
        N=N+1;
    else
        ys=EventWindows(:,E);
        E=E+1;
    end     
%% Getting the Sample differences from the sampled data   
% for j=1:length(ys)-1
% sample_differences(j,runs)= ys(j) - ys(j+1); % taking differecnes of the samples 
% end  
%% Experimenting wihtout doing sample differences 
%Comment this section if you want to use sample differences in the windows
%instead of samples. 
for j=1:length(ys)
sample_differences(j,runs)= ys(j); % just taking the samples 
end  
%% Getting the MLE estimates of the parameter
A = [0,0,-1,1]; %for inequality constraint : 0*Alpha+0*Tau-TauOne+TauTwo<=0
b = 0; 
lb= [-10,0.1,10,1]; % lower bound of the estimates:[alpha tau TauOne TauTwo]
ub= [10,T-.001,100,10]; %upper bound of the estimates
Aeq = []; %for equality constraint
beq = []; %for equality constraint
x=zeros(iter,NumOfParameters);
x0(1,:) = [lb(1),1,lb(3),lb(4)]; % initial estiamtes of the parameters
x0(2,:) = [(ub(1)-lb(1))/3,T/3,(ub(3)-lb(3))/3,(ub(4)-lb(4))/3];
x0(3,:) = [(ub(1)-lb(1))/2,T/2,(ub(3)-lb(3))/2,(ub(4)-lb(4))/2]; 
x0(4,:) = [(ub(1)-lb(1))*3/4,T*3/4,(ub(3)-lb(3))*3/4,(ub(4)-lb(4))*3/4]; 
x0(5,:) = [ub(1),T,ub(3),ub(4)]; 
x0(6,:) = [(ub(1)-lb(1))*5/6,T*5/6,(ub(3)-lb(3))*5/6,(ub(4)-lb(4))*5/6]; 
x0(7,:) = [(ub(1))*2/3,T*2/3,(ub(3))*2/3,(ub(4))*2/3]; 
x0(8,:) = [(ub(1))/4,T/4,(ub(3))/4,(ub(4))/4]; 
x0(9,:) = [(ub(1))/3,T/2,(ub(3))/3,(ub(4))/3]; 
x0(10,:) = [(ub(1))*rand,T/2,(ub(3))*rand,(ub(4))/2]; 
%% When not doing pairwise sample differences in the time windows we pass the derivative as well for optimisation.
for i=1:iter
     if(i>10)
 x0(i,:) = [ub(1)*i/iter,T*i/iter,ub(3)*i/iter,ub(4)*i/iter]; % initial estiamtes of the parameters   
     end
options = optimoptions('fmincon','SpecifyObjectiveGradient',true); 
[x(i,:),val(runs,i),flag(i)] = fmincon(@(t_prime)LikelihoodUnderAlternativeWithDerivative(t_prime,sample_differences(:,runs),noise_var,T),x0(i,:),A,b,Aeq,beq,lb,ub,[],options);
end
%% Doing sample differences in the time windows.
%  for i=1:iter
% % if(i>10)
% % x0(i,:) = [ub(1)*i/iter,T*i/iter,ub(3)*i/iter,ub(4)*i/iter]; % initial estiamtes of the parameters   
% % end
% [x(i,:),val(runs,i),flag(i)] = fmincon(@(t_prime)LikelihoodUnderAlternative(t_prime,sample_differences(:,runs),noise_var,T),x0(i,:),A,b,Aeq,beq,lb,ub);
% end
%% estimates after optimisation.
[~,I]=min(val(runs,:)); % checking for the minm. value of the Likelihood function
mle_estimates(runs,:)=x(I,:); % getting the estimate values(alpha,tau,tauOne,tauTwo) for the corresponding minm. value of the likelihood
alpha_mle(runs)=mle_estimates(runs,1);
tau_mle(runs)=mle_estimates(runs,2);
tauOne_mle(runs)=mle_estimates(runs,3);
tauTwo_mle(runs)=mle_estimates(runs,4);
%% Calculating the Likelihood values under the Null,Alternative hypothesis,deviance while doing pairwise sample differences in the time windows.
%Uncomment this section when using pairwise differenced samples in the
%windows.

% Likelihood_under_Null(runs)= -(length(sample_differences(:,runs))/2)*log(2*pi*(2*noise_var))-(sum((sample_differences(:,runs)).^2))/(2*(2*noise_var));
% Likelihood_under_Alternative(runs) = -val(runs,I);
% %scores(runs)=(Likelihood_under_Alternative(runs) - Likelihood_under_Null(runs)); % Scores are the differences between the Likelihood values
%  scores(runs) = -1*(Likelihood_under_Alternative(runs)); % scores are the negative of the Likelihood values under the alternative hy
%% Calculating the Likelihood values under the Null,Alternative hypothesis,deviance when the samples in the windows are not differenced pairwise. 

 Likelihood_under_Null(runs)= -(length(sample_differences(:,runs))/2)*log(2*pi*(noise_var)) - (sum((sample_differences(:,runs)).^2))/(2*(noise_var));
 Likelihood_under_Alternative(runs) = -val(runs,I);
scores(runs) = -1*(Likelihood_under_Alternative(runs)); % scores are the negative of the Likelihood values under the alternative hypothesis.
 %scores(runs) = (Likelihood_under_Alternative(runs) - Likelihood_under_Null(runs)); % Scores are the differences between the Likelihood values

end
 [pks,loc] = findpeaks(scores); % Threshold is lowered to see if the model actually picks up all the events.
  [X,Y,Threshold,AUC] = perfcurve(labels,scores,1);
%   %% Getting AUC with the peaks 
% NewScores=zeros(1,length(scores));
%   for i=1:length(loc)
%   NewScores(loc(i))=pks(i);
%   end
%% with False positives 
% for i=1:size(X,1) 
%     if(X(i)>=0.30) 
%         I=i; % Index of the point where the FPR is 30% or above
%         break;
%     end
% end
% ClassifierThreshold= Threshold(I); % Threhold score value corresponding to 30% FPR
% L=1;
% for i=1:size(scores,2) 
%     if(scores(i)>ClassifierThreshold)
%         IndexOftheClassifiedEvents(L)=i;
%         L=L+1;
%     end
% end
% loc(:) = IndexOftheClassifiedEvents(:);
 
 end

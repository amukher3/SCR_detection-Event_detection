function[AlphaMLE,TauMLE,TauOneMLE,TauTwoMLE,FirstSample,SclTimeDomain,NumOfWindows,T,SamplingFreqn]=EMForTheProblem(user,DataIndex,Delta);
[Window,FirstSample,SclTimeDomain,T,WindowSize,SamplingFreqn]=MakingWindowsAroundPeaks(user,DataIndex,Delta)
NumIter = 5;  % No. of iterations for the EM
iter=10; % Num of independent initialisations
A=[];b=[];A_prime = [-1,1]; b_prime = [-0.1];
lb= [.01,0.01]; % lower bound of the estimates:[alpha tau]
ub= [5,T-.001]; %upper bound of the estimates
lb_prime= [1,.1]; % lower bound of the estimates:[TauOne TauTwo]
ub_prime= [20,3]; %upper bound of the estimates
Aeq = []; %for equality constraint
beq = []; %for equality constraint
NumOfParameters = 2;
NumOfWindows = size(Window,2);
TauOne(1)=2;
TauTwo(1)=0.75;
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
xPrime=zeros(iter,NumOfParameters);
x0_prime(1,:) = [lb_prime(1)/2,.01]; % initial estiamtes of the parameters
x0_prime(2,:) = [lb_prime(1)/3,1];
x0_prime(3,:) = [lb_prime(1)/5,2];
x0_prime(4,:) = [lb_prime(1)/10,3];
x0_prime(5,:) = [lb_prime(1)/4,4];
x0_prime(6,:) = [ub_prime(1)/6,.01];
x0_prime(7,:) = [ub_prime(1)/ub(1),0.5];
x0_prime(8,:) = [ub_prime(1),1.5];
x0_prime(9,:) = [ub_prime(1),2.5];
x0_prime(10,:) = [ub_prime(1)/2,3.5];
x0_prime(11,:) = [ub_prime(1)/3,4.5];
x0_prime(12,:) = [(lb_prime(1)-ub_prime(1))/2,1.2];
x0_prime(13,:) = [(lb_prime(1)-ub_prime(1))/6,2.2];
x0_prime(14,:) = [(lb_prime(1)- ub_prime(1))/8,3.2];
x0_prime(15,:) = [(lb_prime(1)-ub_prime(1))/3,4.6];
S_d=0;S_r=0;
for Num1=1:NumIter-1
%% E- step :    
for runs=1:NumOfWindows % Num of Windows in a user
  samples(:,runs) = Window(:,runs);
for i=1:iter 
%options = optimoptions('fmincon','SpecifyObjectiveGradient',true); 
[x(i,:),val(runs,i),flag(runs,i)] = fmincon(@(t_prime)LeastSquareEstimates(t_prime,samples(:,runs),T,TauOne(Num1),TauTwo(Num1)),x0(i,:),A,b,Aeq,beq,lb,ub);
end
[fit(runs,Num1),I(Num1)] = min(val(runs,:)); % checking for the minm. value of the Likelihood function
mle_estimates(runs,:) = x(I(Num1),:); % getting the estimate values(alpha,tau,tauOne,tauTwo) for the corresponding minm. value of the likelihood
Alpha(runs,Num1)= mle_estimates(runs,1);
Tau(runs,Num1) = mle_estimates(runs,2);
 end
%% M-step :
samples(:,:) = Window(:,:);
SumSamples(:) = sum(samples(:,:),2);
for i=1:iter 
[x_prime(i,:),val_prime(Num1,i),flag_prime(i)] = fmincon(@(t_prime)LeastSquareEstimates_prime(t_prime,SumSamples(:,:),T,Alpha(:,Num1),Tau(:,Num1)),x0_prime(i,:),A_prime,b_prime,Aeq,beq,lb_prime,ub_prime);
end
[fit_prime(Num1),Iprime(Num1)] = min(val_prime(Num1,:)); % checking for the minm. value of the Likelihood function
mle_estimates_prime(Num1,:)= x_prime(Iprime(Num1),:); % getting the estimate values(alpha,tau,tauOne,tauTwo) for the corresponding minm. value of the likelihood
TauOne(Num1+1) = mle_estimates_prime(Num1,1);
TauTwo(Num1+1) = mle_estimates_prime(Num1,2);
end
AlphaMLE=Alpha(:,Num1);
TauMLE=Tau(:,Num1);
TauOneMLE=TauOne(Num1);
TauTwoMLE=TauTwo(Num1);


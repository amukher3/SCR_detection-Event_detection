% Generative model for undersampled SCR:
% SCR amplitude: alpha ~ Gamma(hyperparam(1),hyperparam(2))
% Bateman decay parameter: tau1 ~ Gamma(hyperparam(3),hyperparam(4))
% Bateman rise parameter: tau2 ~ Gamma(hyperparam(5),hyperparam(6))
% Time shift of SCR: delta ~ Uniform(0,samplePeriod)
% Skin conductance: sc(t) ~ N(Bateman(t),hyperparam(7)) for t =
%   i*samplePeriod, i = 0, 1, ..., nSamples

hyperparam = [2,1,2,2,3,0.2,1e-4];
nSamples = 8;
T = 5;
nWindows = 5000;
globalBateman = false;
% Initial estimates of parameters: [alpha delta tau1 tau2]
paramInit = [2 2.5 3 2];
% Lower and upper bounds of the estimates: [alpha delta tau1 tau2]
lb = [0.01,0.001,0.01,0.01];
ub = [10,T,25,5];

%% Simulate windows
[sc,alpha,delta,tau1,tau2,mu] = simulateWindowsIID(nSamples,T, ...
    hyperparam,nWindows,globalBateman);
% tau1(1)
% tau2(1)
% 
% sseActualParam = sseBateman(sc,T,alpha,delta,tau1,tau2)
% negPostProbActualParam = negPostProbBateman(sc,T,alpha,delta,tau1,tau2, ...
%     hyperparam)

%% Estimation over tau1, tau2 with known alpha, delta
% disp('Optimizing over tau1, tau2')
% 
% % Optimize over global (per-user) tau1, tau2
% [tauMle,sseMle] = fminunc(@(tau) sseBateman(sc,T,alpha,delta,repmat(tau(1), ...
%     nWindows,1),repmat(tau(2),nWindows,1)),paramInit([3 4]));
% tauMle(1)
% tauMle(2)
% sseMle
% 
% [tauMap,negPostProbMap] = fminunc(@(tau) negPostProbBateman(sc,T,alpha,delta, ...
%     repmat(tau(1),nWindows,1),repmat(tau(2),nWindows,1),hyperparam), ...
%     paramInit([3 4]));
% tauMap(1)
% tauMap(2)
% 
% % Try to estimate tau1, tau2 separately for each window
% tauMleAll = zeros(nWindows,2);
% sseMleAll = 0;
% for i = 1:nWindows
%     [tauMleW,sseMleW] = fminunc(@(tau) sseBateman(sc(i,:),T,alpha(i), ...
%         delta(i),tau(1),tau(2)),paramInit([3 4]));
%     tauMleAll(i,:) = tauMleW;
%     sseMleAll = sseMleAll + sseMleW;
% end

%% Try to estimate all 4 parameters separately for each window
disp('Optimizing over all 4 parameters')
% paramMleAll(i,1): estimate for alpha(i)
% paramMleAll(i,2): estimate for delta(i)
% paramMleAll(i,3): estimate for tau1(i)
% paramMleAll(i,4): estimate for tau2(i)
paramMleAll = zeros(nWindows,4);
sseMleAll = zeros(nWindows,1);
exitFlagMle = zeros(nWindows,1);
paramMapAll = zeros(nWindows,4);
negPostProbMapAll = zeros(nWindows,1);
exitFlagMap = zeros(nWindows,1);

parfor i = 1:nWindows
    if mod(i,100) == 0
        disp(['Processing window ' int2str(i)])
    end
    
    options = optimoptions('fmincon','display','off', ...
        'SpecifyObjectiveGradient',false);
    [paramMleW,sseMleAll(i),exitFlagMle(i)] = fmincon(@(param) sseBateman( ...
        sc(i,:),T,param(1),param(2),param(3),param(4)),paramInit, ...
        [],[],[],[],lb,ub,[],options);
    paramMleAll(i,:) = paramMleW;

    options = optimoptions('fminunc','MaxFunctionEvaluations',1000, ...
        'display','off','SpecifyObjectiveGradient',false);
    [paramMapW,negPostProbMapAll(i),exitFlagMap(i)] = fminunc(@(param) ...
        negPostProbBateman(sc(i,:),T,param(1),param(2),param(3),param(4), ...
        hyperparam),paramInit,options);
    paramMapAll(i,:) = paramMapW;
end

%% Summary statistics
% Mean-squared errors
alphaMleMse = mean((alpha-paramMleAll(:,1)).^2);
alphaMapMse = mean((alpha-paramMapAll(:,1)).^2);
deltaMleMse = mean((delta-paramMleAll(:,2)).^2);
deltaMapMse = mean((delta-paramMapAll(:,2)).^2);
tau1MleMse = mean((tau1-paramMleAll(:,3)).^2);
tau1MapMse = mean((tau1-paramMapAll(:,3)).^2);
tau2MleMse = mean((tau2-paramMleAll(:,4)).^2);
tau2MapMse = mean((tau2-paramMapAll(:,4)).^2);

% Standard deviations
alphaMleMseStd = std((alpha-paramMleAll(:,1)).^2)/sqrt(nWindows);
alphaMapMseStd = std((alpha-paramMapAll(:,1)).^2)/sqrt(nWindows);
deltaMleMseStd = std((delta-paramMleAll(:,2)).^2)/sqrt(nWindows);
deltaMapMseStd = std((delta-paramMapAll(:,2)).^2)/sqrt(nWindows);
tau1MleMseStd = std((tau1-paramMleAll(:,3)).^2)/sqrt(nWindows);
tau1MapMseStd = std((tau1-paramMapAll(:,3)).^2)/sqrt(nWindows);
tau2MleMseStd = std((tau2-paramMleAll(:,4)).^2)/sqrt(nWindows);
tau2MapMseStd = std((tau2-paramMapAll(:,4)).^2)/sqrt(nWindows);

fprintf('MSE for each parameter:\n')
fprintf('\\hline\n')
fprintf('Parameter & MLE           & MAP           \\\\\n')
fprintf('\\hline\n')
fprintf('$\\alpha$  & %4.2f \\pm %4.2f & %4.2f \\pm %4.2f \\\\\n', ...
    alphaMleMse,alphaMleMseStd,alphaMapMse,alphaMapMseStd)
fprintf('$\\delta$  & %4.2f \\pm %4.2f & %4.2f \\pm %4.2f \\\\\n', ...
    deltaMleMse,deltaMleMseStd,deltaMapMse,deltaMapMseStd)
fprintf('\\tau_1    & %4.2f \\pm %4.2f & %4.2f \\pm %4.2f \\\\\n', ...
    tau1MleMse,tau1MleMseStd,tau1MapMse,tau1MapMseStd)
fprintf('\\tau_2    & %4.2f \\pm %4.2f & %4.2f \\pm %4.2f \\\\\n', ...
    tau2MleMse,tau2MleMseStd,tau2MapMse,tau2MapMseStd)
fprintf('\\hline\n')

%% Scatter plots of parameter estimates and actual values
figure(1)
subplot(221), scatter(alpha,paramMleAll(:,1),'.')
axis equal
xlabel alpha
ylabel MLE
subplot(222), scatter(delta,paramMleAll(:,2),'.')
axis equal
xlabel delta
ylabel MLE
subplot(223), scatter(tau1,paramMleAll(:,3),'.')
axis equal
xlabel tau1
ylabel MLE
subplot(224), scatter(tau2,paramMleAll(:,4),'.')
axis equal
xlabel tau2
ylabel MLE

figure(2)
subplot(221), scatter(alpha,paramMapAll(:,1),'.')
axis equal
xlabel alpha
ylabel MAP
subplot(222), scatter(delta,paramMapAll(:,2),'.')
axis equal
xlabel delta
ylabel MAP
subplot(223), scatter(tau1,paramMapAll(:,3),'.')
axis equal
xlabel tau1
ylabel MAP
subplot(224), scatter(tau2,paramMapAll(:,4),'.')
axis equal
xlabel tau2
ylabel MAP

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
nWindows = 5;
globalBateman = false;
% Initial estimates of parameters: [alpha delta tau1 tau2]
paramInit = [2 2.5 3 2];
% Lower and upper bounds of the estimates: [alpha delta tau1 tau2]
lb = [0.01,0.001,0.01,0.01];
ub = [10,T,25,5];

%% Simulate windows
[sc,alpha,delta,tau1,tau2,mu] = simulateWindowsIID(nSamples,T, ...
    hyperparam,nWindows,globalBateman);
tau1(1)
tau2(1)

sseActualParam = sseBateman(sc,T,alpha,delta,tau1,tau2)
negPostProbActualParam = negPostProbBateman(sc,T,alpha,delta,tau1,tau2, ...
    hyperparam)

%% Grid search over tau1, tau2
disp('Optimizing over tau1, tau2 by grid search')
% tau1Range = 0.1:0.1:10;
% tau2Range = 0.1:0.1:5;
% sse = zeros(length(tau1Range),length(tau2Range));
% negPostProb = zeros(length(tau1Range),length(tau2Range));
% for i = 1:length(tau1Range)
%     if mod(tau1Range(i),1) == 0
%         disp(['Grid search: tau1 = ' num2str(tau1Range(i))])
%     end
%     
%     for j = 1:length(tau2Range)
%         sse(i,j) = sseBateman(sc,T,alpha,delta,tau1Range(i)*ones(nWindows, ...
%             1),tau2Range(j)*ones(nWindows,1));
%         negPostProb(i,j) = negPostProbBateman(sc,T,alpha,delta,tau1Range(i) ...
%             *ones(nWindows,1),tau2Range(j)*ones(nWindows,1),hyperparam);
%     end
% end
% [sseMin,indexMin] = min(sse(:));
% [iMin,jMin] = ind2sub([length(tau1Range),length(tau2Range)],indexMin);
% tau1Range(iMin)
% tau2Range(jMin)
% sseMin
% 
% [nppMin,indexMin] = min(negPostProb(:));
% [iMin,jMin] = ind2sub([length(tau1Range),length(tau2Range)],indexMin);
% tau1Range(iMin)
% tau2Range(jMin)
% nppMin

% Grid search over all 4 parameters
alphaRange = logspace(-1,1,20);
deltaRange = 0:0.1:T;
tau1Range = logspace(-1,log10(25),20);
tau2Range = logspace(-1,log10(5),20);
sse = zeros(length(alphaRange),length(deltaRange),length(tau1Range), ...
    length(tau2Range));
negPostProb = zeros(length(alphaRange),length(deltaRange),length(tau1Range), ...
    length(tau2Range));
window = 1;
for i1 = 1:length(alphaRange)
    disp(['Grid search: alpha = ' num2str(alphaRange(i1))])
    for i2 = 1:length(deltaRange)
        for i3 = 1:length(tau1Range)
            for i4 = 1:length(tau2Range)
                sse(i1,i2,i3,i4) = sseBateman(sc(window,:),T,alphaRange(i1), ...
                    deltaRange(i2),tau1Range(i3),tau2Range(i4));
                negPostProb(i1,i2,i3,i4) = negPostProbBateman(sc(window,:), ...
                    T,alphaRange(i1),deltaRange(i2),tau1Range(i3), ...
                    tau2Range(i4),hyperparam);
            end
        end
    end
end
[sseMin,indexMin] = min(sse(:));
[i1Min,i2Min,i3Min,i4Min] = ind2sub([length(alphaRange),length(deltaRange), ...
    length(tau1Range),length(tau2Range)],indexMin);
alphaRange(i1Min)
deltaRange(i2Min)
tau1Range(i3Min)
tau2Range(i4Min)
sseMin

[nppMin,indexMinP] = min(negPostProb(:));
[i1MinP,i2MinP,i3MinP,i4MinP] = ind2sub([length(alphaRange), ...
    length(deltaRange),length(tau1Range),length(tau2Range)],indexMinP);
alphaRange(i1MinP)
deltaRange(i2MinP)
tau1Range(i3MinP)
tau2Range(i4MinP)
nppMin


%% Plot grid search result
% figure(1)
% imagesc(sse)
% colorbar
% ylabel tau1
% xlabel tau2
% title(['SSE: actual tau1=' num2str(tau1(1)) ', tau2=' num2str(tau2(2))])
% set(gca,'yticklabel',tau1Range(get(gca,'ytick')))
% set(gca,'xticklabel',tau1Range(get(gca,'xtick')))
% 
% figure(2)
% imagesc(negPostProb)
% colorbar
% ylabel tau1
% xlabel tau2
% title(['-log(posterior probability): actual tau1=' num2str(tau1(1)) ...
%     ', tau2=' num2str(tau2(2))])
% set(gca,'yticklabel',tau1Range(get(gca,'ytick')))
% set(gca,'xticklabel',tau1Range(get(gca,'xtick')))

% Plots of SSE
figure
subplot(231)
h = pcolor(alphaRange,deltaRange,log(squeeze(sse(:,:,i3Min,i4Min)))');
h.EdgeColor = 'none';
colorbar
title(['Actual alpha=' num2str(alpha(window)) ', delta=' ...
    num2str(delta(window))])
xlabel alpha
ylabel delta

subplot(232)
h = pcolor(alphaRange,tau1Range,log(squeeze(sse(:,i2Min,:,i4Min)))');
h.EdgeColor = 'none';
colorbar
title(['Actual alpha=' num2str(alpha(window)) ', tau1=' ...
    num2str(tau1(window))])
xlabel alpha
ylabel tau1

subplot(233)
h = pcolor(alphaRange,tau2Range,log(squeeze(sse(:,i2Min,i3Min,:)))');
h.EdgeColor = 'none';
colorbar
title(['Actual alpha=' num2str(alpha(window)) ', tau2=' ...
    num2str(tau2(window))])
xlabel alpha
ylabel tau2

subplot(234)
h = pcolor(deltaRange,tau1Range,log(squeeze(sse(i1Min,:,:,i4Min)))');
h.EdgeColor = 'none';
colorbar
title(['Actual delta=' num2str(delta(window)) ', tau1=' ...
    num2str(tau1(window))])
xlabel delta
ylabel tau1

subplot(235)
h = pcolor(deltaRange,tau2Range,log(squeeze(sse(i1Min,:,i3Min,:)))');
h.EdgeColor = 'none';
colorbar
title(['Actual delta=' num2str(delta(window)) ', tau2=' ...
    num2str(tau2(window))])
xlabel delta
ylabel tau2

subplot(236)
h = pcolor(tau1Range,tau2Range,log(squeeze(sse(i1Min,i2Min,:,:)))');
h.EdgeColor = 'none';
colorbar
title(['Actual tau1=' num2str(tau1(window)) ', tau2=' ...
    num2str(tau2(window))])
xlabel tau1
ylabel tau2

% Plots of negative log-posterior
figure
subplot(231)
h = pcolor(alphaRange,deltaRange,log(squeeze(sse(:,:,i3MinP,i4MinP)))');
h.EdgeColor = 'none';
colorbar
title(['Actual alpha=' num2str(alpha(window)) ', delta=' ...
    num2str(delta(window))])
xlabel alpha
ylabel delta

subplot(232)
h = pcolor(alphaRange,tau1Range,log(squeeze(negPostProb(:,i2MinP,:,i4MinP)))');
h.EdgeColor = 'none';
colorbar
title(['Actual alpha=' num2str(alpha(window)) ', tau1=' ...
    num2str(tau1(window))])
xlabel alpha
ylabel tau1

subplot(233)
h = pcolor(alphaRange,tau2Range,log(squeeze(negPostProb(:,i2MinP,i3MinP,:)))');
h.EdgeColor = 'none';
colorbar
title(['Actual alpha=' num2str(alpha(window)) ', tau2=' ...
    num2str(tau2(window))])
xlabel alpha
ylabel tau2

subplot(234)
h = pcolor(deltaRange,tau1Range,log(squeeze(negPostProb(i1MinP,:,:,i4MinP)))');
h.EdgeColor = 'none';
colorbar
title(['Actual delta=' num2str(delta(window)) ', tau1=' ...
    num2str(tau1(window))])
xlabel delta
ylabel tau1

subplot(235)
h = pcolor(deltaRange,tau2Range,log(squeeze(negPostProb(i1MinP,:,i3MinP,:)))');
h.EdgeColor = 'none';
colorbar
title(['Actual delta=' num2str(delta(window)) ', tau2=' ...
    num2str(tau2(window))])
xlabel delta
ylabel tau2

subplot(236)
h = pcolor(tau1Range,tau2Range,log(squeeze(negPostProb(i1MinP,i2MinP,:,:)))');
h.EdgeColor = 'none';
colorbar
title(['Actual tau1=' num2str(tau1(window)) ', tau2=' ...
    num2str(tau2(window))])
xlabel tau1
ylabel tau2


function [sc,alpha,delta,tau1,tau2,mu] = simulateWindowsIID(nSamples...
    samplePeriod,hyperparam,nWindows,globalBateman)

% Generative model for undersampled SCR:
% SCR amplitude: alpha ~ Gamma(hyperparam(1),hyperparam(2))
% Bateman decay parameter: tau1 ~ Gamma(hyperparam(3),hyperparam(4))
% Bateman rise parameter: tau2 ~ Gamma(hyperparam(5),hyperparam(6))
% Time shift of SCR: delta ~ Uniform(0,samplePeriod)
% Skin conductance: sc(t) ~ N(Bateman(t),hyperparam(7)) for t =
%   i*samplePeriod, i = 0, 1, ..., nSamples

if nargin < 5
    globalBateman = true;
end

% Times at which samples are drawn in each window
t = 0:samplePeriod:samplePeriod*(nSamples-1);

% Bateman bi-exponential functional form to model SCR
% Sample parameter values

delta = unifrnd(0,samplePeriod,nWindows,1);
alpha = gamrnd(hyperparam(1),hyperparam(2),nWindows,1);

% Generate a single set of tau1,tau2 parameters for the Bateman function if
% globalBateman is true

if globalBateman == true
    tau1 = gamrnd(hyperparam(3),hyperparam(4))*ones(nWindows,1);
    tau2 = gamrnd(hyperparam(5),hyperparam(6))*ones(nWindows,1);
else
    tau1 = gamrnd(hyperparam(3),hyperparam(4),nWindows,1);
    tau2 = gamrnd(hyperparam(5),hyperparam(6),nWindows,1);
end

% Sample skin conductance windows
mu = zeros(nWindows,nSamples);  % Noiseless samples

for i = 1:nWindows
    mu(i,:) = bateman(t,alpha(i),delta(i),tau1(i),tau2(i));
end
%sc = mu + sqrt(hyperparam(7))*randn(nWindows,nSamples);
sc = mu + sqrt(hyperparam(7));
end


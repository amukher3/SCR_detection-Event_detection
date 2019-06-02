function [negLogProb,gradient] = negPostProbBateman(sc,samplePeriod,alpha, ...
    delta,tau1,tau2,hyperparam)
% Generative model for undersampled SCR:
% SCR amplitude: alpha ~ Gamma(hyperparam(1),hyperparam(2))
% Bateman decay parameter: tau1 ~ Gamma(hyperparam(3),hyperparam(4))
% Bateman rise parameter: tau2 ~ Gamma(hyperparam(5),hyperparam(6))
% Time shift of SCR: delta ~ Uniform(0,samplePeriod)
% Skin conductance: sc(t) ~ N(Bateman(t),hyperparam(7)) for t =
%   i*samplePeriod, i = 0, 1, ..., nSamples

[nWindows,nSamples] = size(sc);

sse = sseBateman(sc,samplePeriod,alpha,delta,tau1,tau2);
logPriorAlpha = log(gampdf(alpha,hyperparam(1),hyperparam(2)));
logPriorDelta = log(unifpdf(delta,0,samplePeriod));
logPriorTau1 = log(gampdf(tau1,hyperparam(3),hyperparam(4)));
logPriorTau2 = log(gampdf(tau2,hyperparam(5),hyperparam(6)));

logProb = -sse/(2*hyperparam(7)) + sum(logPriorAlpha + logPriorDelta + ...
    logPriorTau1 + logPriorTau2);
negLogProb = -logProb;

% Compute gradient as second output if requested
if nargout > 1    
    gradient = zeros(4,1);
    for i = 1:nSamples
        t = (i-1)*samplePeriod;
        dsse_dmu = -2*(sc(i) - bateman(t,alpha,delta,tau1,tau2));
        dmu_dalpha = bateman(t,alpha,delta,tau1,tau2) / alpha;
        if t >= delta
            dmu_ddelta = alpha*(exp(-(t-delta)/tau1)/tau1 ...
                - exp(-(t-delta)/tau1)/tau2);
            dmu_dtau1 = alpha/(tau1^2) * (t-delta) * exp(-(t-delta)/tau1);
            dmu_dtau2 = -alpha/(tau2^2) * (t-delta) * exp(-(t-delta)/tau2);
        else
            dmu_ddelta = 0;
            dmu_dtau1 = 0;
            dmu_dtau2 = 0;
        end
        dsse_dalpha = dsse_dmu * dmu_dalpha;
        dsse_ddelta = dsse_dmu * dmu_ddelta;
        dsse_dtau1 = dsse_dmu * dmu_dtau1;
        dsse_dtau2 = dsse_dmu * dmu_dtau2;
        
        dlogProb_dalpha = -dsse_dalpha/(2*hyperparam(7)) + (hyperparam(1) ...
            - 1) / alpha - 1/hyperparam(2);
        dlogProb_ddelta = -dsse_ddelta/(2*hyperparam(7));
        dlogProb_dtau1 = -dsse_dtau1/(2*hyperparam(7)) + (hyperparam(3) ...
            - 1) / tau1 - 1/hyperparam(4);
        dlogProb_dtau2 = -dsse_dtau2/(2*hyperparam(7)) + (hyperparam(5) ...
            - 1) / tau2 - 1/hyperparam(6);
        
        % Subtract instead of add because the gradient computed is for the
        % log posterior probability, not the negative
        gradient = gradient - [dlogProb_dalpha;
            dlogProb_ddelta;
            dlogProb_dtau1;
            dlogProb_dtau2];
    end
end

end


function [sse,gradient] = sseBateman(sc,samplePeriod,alpha,delta,tau1,tau2)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[nWindows,nSamples] = size(sc);
tVec = 0:samplePeriod:samplePeriod*(nSamples-1);

sse = 0;
for i = 1:nWindows
    mu = bateman(tVec,alpha(i),delta(i),tau1(i),tau2(i));
    sse = sse + sum((sc(i,:) - mu).^2);
end

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
        gradient = gradient + [dsse_dmu * dmu_dalpha;
            dsse_dmu * dmu_ddelta;
            dsse_dmu * dmu_dtau1;
            dsse_dmu * dmu_dtau2];
    end
end

end


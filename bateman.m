function output = bateman(t,alpha,delta,tau1,tau2)
%bateman Bateman bi-exponential form for skin conductance response (SCR)
%   bateman(t,alpha,delta,tau1,tau2) returns the value of the shifted and
%   scaled Bateman bi-exponential form for an SCR
%       b(t) = alpha(e^-(t-delta)/tau1 - e^-(t-delta)/tau2) for t >= delta
%       b(t) = 0                                            for t < delta
%   at the vector of times given by t

[nWindows,nSamples] = size(t);
validTimes = (t >= delta);
output = zeros(nWindows,nSamples);
if ~isempty(t(validTimes))
    output(validTimes) = alpha.*(exp(-(t(validTimes)-delta)./tau1) ...
        - exp(-(t(validTimes)-delta)./tau2));
end

end


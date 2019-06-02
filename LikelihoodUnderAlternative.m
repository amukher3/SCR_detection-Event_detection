function L = LikelihoodUnderAlternative(x,sample_differences,noise_var,T)
 Lhood=0;
for h=1:length(sample_differences)
%% with pariwise differences of the samples
% mu_1(h)= x(1)*(exp((x(2)-(h+1)*T)/x(3)) - exp((x(2)-(h+1)*T)/x(4))).*heaviside(-x(2)+(h+1)*T);
% mu_2(h) = x(1)*(exp((x(2)-(h)*T)/x(3)) - exp((x(2)-(h)*T)/x(4))).*heaviside(-x(2)+(h)*T);
% mu(h)=mu_1(h) - mu_2(h);
% Lhood = Lhood + (-1/2*log(2*pi*2*(noise_var)) - (((sample_differences(h)-mu(h))^2)/(2*(2*noise_var))));

%% no pairwise difference of the samples
mu(h)= x(1)*(exp((x(2)-(h)*T)/x(3))- exp((x(2)-(h)*T)/x(4)));
Lhood = Lhood + (-1/2*log(2*pi*(noise_var)) -(((sample_differences(h)-mu(h))^2)/(2*(noise_var))));
end
L=-(Lhood);
end





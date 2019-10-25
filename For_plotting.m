close all
SamplingFreq=32;
T=5;
user=user3;
DownSamFact=T*SamplingFreq;
limit=2500;
Delta=1;
SampleInstants=Delta*SamplingFreq:DownSamFact:length(user);
DownSamSig=user(SampleInstants);
plot(user)
hold on
plot(SampleInstants,DownSamSig,'o')
xlabel('Samples')
ylabel('Skin Conductance')
xlim([1000 limit])
% subplot(2,2,1)
% plot(user)
% hold on
% plot(SampleInstants,DownSamSig,'o')
% xlabel('Samples')
% ylabel('Skin Conductance')
% xlim([1 limit])
% Delta=2;
% SampleInstants=Delta*SamplingFreq:DownSamFact:length(user);
% DownSamSig=user(SampleInstants);
% subplot(2,2,2)
% plot(user)
% hold on
% plot(SampleInstants,DownSamSig,'o')
% xlim([1 limit])
% xlabel('Samples')
% ylabel('Skin Conductance')
% Delta=3;
% SampleInstants=Delta*SamplingFreq:DownSamFact:length(user);
% DownSamSig=user(SampleInstants);
% subplot(2,2,3)
% plot(user)
% hold on
% plot(SampleInstants,DownSamSig,'o')
% xlim([1 limit])
% xlabel('Samples')
% ylabel('Skin Conductance')
% Delta=4;
% SampleInstants=Delta*SamplingFreq:DownSamFact:length(user);
% DownSamSig=user(SampleInstants);
% subplot(2,2,4)
% plot(user)
% hold on
% plot(SampleInstants,DownSamSig,'o')
% xlim([1 limit])
% xlabel('Samples')
% ylabel('Skin Conductance')
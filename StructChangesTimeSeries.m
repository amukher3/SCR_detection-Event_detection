%YogaTimeSeries=[YogaSeries_RR(:,2);YogaSeries_RR3(:,2);YogaSeries_RR4(:,2)];
findchangepts(TaiChiSeries,'MaxNumChanges',7,'Statistic','mean');
% YogaLen=length(YogaTimeSeries);
% t=1:1:YogaLen;
% SubInit=[1,1585,3215];
% StimIns=[570,2305,3659];
% unitstep=t>=SubInit(:) & t<=SubInit(:)+10;
% hold on
% h1=plot(t,140*unitstep,'--');
% set(h1,'color','red')
% unitstep=t>=StimIns(:) & t<=StimIns(:)+10;
% h2=plot(t,140*unitstep);
% set(h2,'color','green')
% legend('HR trace','New Subject','Stimulus Instants')
% xlabel('Samples');
% ylabel('Heart rate(BPM)');
% title('HR time series for Yoga');
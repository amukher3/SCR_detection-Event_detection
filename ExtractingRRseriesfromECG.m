clear all
close all

NumUsers=10;
for UserIdx=1:NumUsers
    
if (UserIdx<=9)
        formatSpec =...
            'C:\\Users\\Abhishek Mukherjee\\Downloads\\f2y0%dm';
        Num=UserIdx;
        str=sprintf(formatSpec,Num);
        load(str);
        
%         %Using RDMAT to extract the signal
%         [tm,signal,Fs,labels]=rdmat(str);
    else
        
        formatSpec =...
            'C:\\Users\\Abhishek Mukherjee\\Downloads\\f2y%dm';
        Num=UserIdx;
        str=sprintf(formatSpec,Num);
        load(str);
        
%         %Using RDMAT to extract the signal
%         [tm,signal,Fs,labels]=rdmat(str);
        
end
  

    
%%% Extracting the ECG & respiratory data %%%
 ppgData=val(2,:);
 respData{UserIdx}=val(1,:);
 
[~,Idx{UserIdx}]=findpeaks(ppgData,'MinPeakDistance',100,...
'MinPeakProminence',0.5);

%%% Creating the TimeDiff Series %%

  for j=1:size(Idx,2)
      
 temp=Idx{j};
 
 for i=1:length(temp)-1
     
     % samples difference between peaks
     DiffTemp(i)=temp(i+1)-temp(i); 
     
     xIdx(i)=temp(i);
     yIdx(temp(i))=DiffTemp(i);
     
 end
 diffVect{j}=DiffTemp;
 
 % Sampling frequency for the PPG signal
 Fs=250; % Sampling frequency of the ECG signal
 
 %InterpolatingFactor for Sinc Interpolation
 IntrpFact=round(length(ppgData)/length(DiffTemp)); 
 
%%% Uncomment for Linear Interpolation
%  xv=1:1:length(DiffTemp);
%  xq=1:1:length(Data.RR(1,:));
%  RR_series{j} = interp1(xv,DiffTemp,xq);

% Interpolating the extracted RR series
% Sinc interpolation..
 RR_series{j} = interp(DiffTemp/Fs,IntrpFact);

 clear temp diffTemp;
 
  end

TempPrime=cell2mat(RR_series(UserIdx));

% figure;
% plot(TempPrime)

[~,IdxPrime]=find(TempPrime<0.3);
[~,IdxDoublePrime]=find(TempPrime>1.4);

%%% Percentage deviation in the expected values 
%%% of the extracted RR series 
%%% Expected IBI range 0.3-1.4 secs.

PercentFac(UserIdx)= (length(IdxPrime)+length(IdxDoublePrime))/...
    length(TempPrime)*100;


 end
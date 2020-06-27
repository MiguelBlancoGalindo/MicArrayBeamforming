function R_xx=CalcRxx(Audio,Config,ArrayMan)

nfft=Config.Filter.Nfft;
nFreqs=length(ArrayMan.f);
nMics=length(Config.Array.MicPos);

% calculate the expectation of the received signals
maxsamples=500;
SegSamples=Config.Filter.SegSamples;
DataSize=size(Audio);
StartSamples=1:SegSamples:max(DataSize);
if length(StartSamples)>maxsamples
   StartSamples(maxsamples+1:end)=[]; 
end
disp(['first ' num2str(length(StartSamples)*SegSamples/16000) 's used']);
R_x=zeros(nMics,nMics,length(StartSamples)-1,nFreqs);

if strcmp(Config.Filter.FreqMode,'discrete')
    pickfftoutput=1;
    fftfreqs=((0:nfft/2-1)./nfft)*Config.Filter.Fs;
    freqinds=zeros(1,length(ArrayMan.f));
    for kk=1:length(ArrayMan.f)
        [~,freqinds(kk)]=min(abs(ArrayMan.f(kk)-fftfreqs));
    end
else
    pickfftoutput=0;
end

if size(Audio,1)>size(Audio,2)
    Audio=Audio.';
end
for j=2:length(StartSamples)
    AudioSeg=Audio(:,StartSamples(j-1):StartSamples(j)-1);
    
    DataFD=fft(AudioSeg,nfft,2);
    
    if pickfftoutput
        DataFD=DataFD(:,freqinds);
    end
        
    for k=1:nFreqs
        R_x(:,:,j-1,k)=DataFD(:,k)*DataFD(:,k)';
    end
end
%toc
R_xx=squeeze(mean(R_x,3));
disp('calculated the Rxx matrix estimate');


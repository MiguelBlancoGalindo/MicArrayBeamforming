function x = SynthesiseMicSignals(Config, ArrayMan)
% function x = SynthesiseMicSignals(Config, ArrayMan)
% function that synthesises the microphone signals captured by a virtual
% array from the settings in Config.
%
% input arguments:
%   Config: configuration structure containing all the settings including
%       paths to input signals, positions of those signals in space and 
%       type of array.
%   ArrayMan: array manifold structure containing the transfer functions to
%       synthesise the signals at the array. 
%
% output arguments:
%   x: microphone array signals.

if size(Config.Scene.Targets.FilePath,1) ~= length(Config.Scene.Targets)
    error('number of targets different from number of specified audio file paths');
elseif size(Config.Scene.Interferers.FilePath,1) ~= length(Config.Scene.Interferers)
    error('number of interferers different from number of specified audio file paths');
else
    TarAngles = Config.Scene.Targets.Angle;
    IntAngles = Config.Scene.Interferers.Angle;
    nTars = size(Config.Scene.Targets.Angle,1);
    nInts = size(Config.Scene.Interferers.Angle,2)/2;
    if nTars~=length(Config.Scene.Targets.Gain)
        error('number of target angles not matching its number of gains');
    elseif nInts~=length(Config.Scene.Interferers.Gain)
        error('number of interferer angles not matching its number of gains');
    end
    
    if ~isfield(Config.Scene.Targets,'AngleInd')
    %find indices of target and interferers
    for iTar=1:nTars
        [~,~,Config.Scene.Targets.AngleInd(iTar)] = getNearestAngle(Config,TarAngles(iTar,1),TarAngles(iTar,2));
    end
    end
    if ~isfield(Config.Scene.Interferers,'AngleInd')
    for iTar=1:nTars
        for iInt=1:nInts            
            [~,~,Config.Scene.Interferers.AngleInd(iTar,iInt)] = getNearestAngle(Config,IntAngles(iTar,iInt*2-1),IntAngles(iTar,iInt*2));
        end
    end
    end
    nSimErrors = length(ArrayMan.Monitor);
    if isfield(Config,'Temp') && ~isempty(Config.Temp) && isfield(Config.Temp,'Sims') && ~isempty(Config.Temp.Sims)
        doSim = zeros(length(ArrayMan.Monitor),1);
        doSim(Config.Temp.Sims) = ones(length(Config.Temp.Sims),1);
    else
        doSim = ones(length(ArrayMan.Monitor),1);
    end
    nSims = sum(doSim==1);
    
    sTar = cell(nTars,1);
    sInt = cell(nInts,1);
    if iscell(Config.Scene.Targets.FilePath)
        tarIsCell = true;
    elseif ischar(Config.Scene.Targets.FilePath)
        tarIsCell = false;
    end
    if iscell(Config.Scene.Interferers.FilePath)
        intIsCell = true;
    elseif ischar(Config.Scene.Interferers.FilePath)
        intIsCell = false;
    end
    
    for iTar=1:nTars
        if ~isfield(Config.Scene.Targets,'SampleRange')
            if tarIsCell;  [sTar{iTar},fs] = audioread(Config.Scene.Targets.FilePath{iTar});
            else [sTar{iTar},fs] = audioread(Config.Scene.Targets.FilePath); end
        else
            if tarIsCell;  [sTar{iTar},fs] = audioread(Config.Scene.Targets.FilePath{iTar}, Config.Scene.Targets.Segment{iTar});
            else [sTar{iTar},fs] = audioread(Config.Scene.Targets.FilePath, Config.Scene.Targets.SampleRange); end
        end
        rows = size(sTar{iTar},1); columns = size(sTar{iTar},2);
        if rows<columns; sTar{iTar}=sTar{iTar}.'; end            
        sTar{iTar} = sTar{iTar}(:,1);
        if fs~=Config.Filter.Fs; sTar{iTar} = resample(sTar{iTar},Config.Filter.Fs,fs); end
        Nconv = length(sTar{1});
    end    
    
    for iInt=1:nInts
        if ~isfield(Config.Scene.Interferers,'SampleRange')
            if intIsCell;  [sInt{iInt},fs] = audioread(Config.Scene.Interferers.FilePath{iInt});
            else [sInt{iInt},fs] = audioread(Config.Scene.Interferers.FilePath); end
        else
            if intIsCell;  [sInt{iInt},fs] = audioread(Config.Scene.Interferers.FilePath{iInt}, Config.Scene.Interferers.SampleRange{iInt});
            else [sInt{iInt},fs] = audioread(Config.Scene.Interferers.FilePath, Config.Scene.Interferers.SampleRange); end
        end
        rows = size(sInt{iInt},1); columns = size(sInt{iInt},2);
        if rows<columns; sInt{iInt}=sInt{iInt}.'; end            
        sInt{iInt} = sInt{iInt}(:,1);
        if fs~=Config.Filter.Fs; sInt{iInt} = resample(sInt{iInt},Config.Filter.Fs,fs); end
    end
    x = zeros(Nconv,Config.Array.M,nSims,nTars);

    if isfield(Config.Scene,'Diffuse') && isfield(Config.Scene.Diffuse.Coeff,'Diffuse') && ...
        Config.Scene.Diffuse.Coeff~=0 && ~isempty(Config.Scene.Diffuse.FilePath)
        doDiffuse = true;
        sDiff = cell(nTars,1);
        xDiff = zeros(Config.Array.M,Config.Filter.Nfft,nSimErrors,nTars);
        L = size(ArrayMan.Monitor(1).a,2).*size(ArrayMan.Monitor(1).a,3);
    else 
        doDiffuse = false;
        xDiff = 0;
    end
    
    is=0;
    for iSim=1:nSimErrors 
        if doSim(iSim)
        is = is+1;
        if doDiffuse
            Adiff=sum(reshape(ArrayMan.Monitor(iSim).a,Config.Array.M,L,Config.Filter.Nfft),2)/L;
            AdiffTD=GetTDArrayMan(Adiff,Config.Array.M,Config);
        end
        for iTar=1:nTars
            aTar = squeeze(ArrayMan.Monitor(iSim).a(:,Config.Scene.Targets.AngleInd(iTar),:));
            if size(aTar,1)>size(aTar,2); aTar = aTar.'; end
            aTarTD=GetTDArrayMan(aTar,Config.Array.M,Config);
            xTar = fftfilt(aTarTD,Config.Scene.Targets.Gain(iTar).*repmat(sTar{iTar},1,Config.Array.M));
            Nconv = length(sInt{1});
            xInt = zeros(Nconv,Config.Array.M,nInts);
            for iInt=1:nInts               
                aInt = squeeze(ArrayMan.Monitor(iSim).a(:,Config.Scene.Interferers.AngleInd(iTar,iInt),:));
                if size(aInt,1)>size(aInt,2); aInt = aInt.'; end
                aIntTD=GetTDArrayMan(aInt,Config.Array.M,Config);
                xInt(:,:,iInt) = fftfilt(aIntTD,Config.Scene.Interferers.Gain(iTar,iInt).*repmat(sInt{iInt},1,Config.Array.M));
            end
            xInt = sum(xInt,3);
        
            if doDiffuse
                sDiff{iTar} = audioread(Config.Scene.Diffuse.FilePath{iTar}).';
                diffGain = sqrt(Config.Scene.Diffuse.Coeff(iTar).*(Config.Scene.Targets.Gain(iTar).^2+sum(Config.Scene.Targets.Gain(iTar,:).^2)));
                xDiff = fftfilt(AdiffTD,diffGain.*sDiff{iTar});
            end
            x1 = xTar + xInt + xDiff;
            x(:,:,is,iTar) = 0.999.*x1./max(max(abs(x1)));
        end
        end
    end   
    
end
    
    
end
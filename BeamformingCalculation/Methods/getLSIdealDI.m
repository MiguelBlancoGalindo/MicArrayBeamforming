function DITarget = getLSIdealDI(Config)
% function DITarget = getLSIdealDI(Config)
% function that obtains the directivity index of the target directivity
% function for the least-squares beamformer based on the settings specified
% in Config.
%
% input argument:
%   Config: configuration struct containing all the settings.
%
% output argument:
%   d: target directivity function.

if ~strcmpi(Config.Filter.Method,'ls')
    error('Unable to compute ideal beampattern as beamforming method is not least squares');
else
    dd = getLSIdealPattern(Config);
    targetAz = Config.Filter.CurrTarget(2);
    targetEl = Config.Filter.CurrTarget(1);
    [~,~,targetInd] = getNearestAngle(Config,targetEl,targetAz);

    DITarget = 10*log10(abs(dd(targetInd)).^2./mean(abs(dd).^2));

end

end
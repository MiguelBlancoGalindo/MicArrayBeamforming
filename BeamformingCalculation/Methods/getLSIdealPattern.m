function d = getLSIdealPattern(Config)
% function d = getTargetPattern(Config)
% function that obtains the target directivity function from the settings
% specified in Config. It uses getTargetPattern(Config) and it is not
% deleted for legacy purposes to be able to run scripts that relied on this
% function.
%
% input argument:
%   Config: configuration struct containing all the settings.
%
% output argument:
%   d: target directivity function.
d = getTargetPattern(Config);

end
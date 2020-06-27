function Config=checkConfigFilters(Config)
% function Config=checkConfigFilters(Config)
% function that checks that the configuration parameters of the beamforming
% calculation are not missing or incorrect. 
%
% input and output arguments: 
%   Config: configuration structure containing all the settings.

LSNMethods=0;
for imethod=1:length(Config.Filter.Method)
    bftype=Config.Filter.Method{imethod};
    if strcmp(bftype,'ls') || strcmp(bftype,'mn')
        if ischar(Config.Filter.(sprintf(bftype)).mode); Config.Filter.(sprintf(bftype)).mode = {Config.Filter.(sprintf(bftype)).mode}; end
        LSNMethods = LSNMethods+1;
    end
end
if isfield(Config.Filter,'ls') 
    LSNModes = length(Config.Filter.ls.mode);
    LSNOrders = length(Config.Filter.ls.order);
    if LSNModes~=LSNMethods && LSNModes==1
        Config.Filter.ls.mode = repmat(Config.Filter.ls.mode,LSNMethods,1);
    end
    if LSNOrders~=LSNMethods && LSNOrders==1
        Config.Filter.ls.order = repmat(Config.Filter.ls.order,LSNMethods,1);
    end
    if LSNOrders~=LSNModes && LSNModes==1
        Config.Filter.ls.mode = repmat(Config.Filter.ls.mode,LSNOrders,1);
    end
end
if isfield(Config.Filter,'mn')
    LSNModes = length(Config.Filter.mn.mode);
    LSNOrders = length(Config.Filter.mn.order);
    if LSNModes~=LSNMethods && LSNModes==1
        Config.Filter.mn.mode = repmat(Config.Filter.mn.mode,LSNMethods,1);
    end
    if LSNOrders~=LSNMethods && LSNOrders==1
        Config.Filter.mn.order = repmat(Config.Filter.mn.order,LSNMethods,1);
    end
    if LSNOrders~=LSNModes && LSNModes==1
        Config.Filter.ls.mode = repmat(Config.Filter.ls.mode,LSNOrders,1);
    end
end


end

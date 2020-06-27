function [BFWeights,EvalStruct,ArrayMan,Rxx]=RunBFTests(Config,ArrayMan,SetupAudio,Rxx)
% [BFWeights,EvalStruct,ArrayMan,Rxx]=RunBFTests(Config,ArrayMan,SetupAudio,Rxx)
%
% Run and evaluate a set of beamformer tests based on the input Config
% struct
%
% ArrayMan can be given as input calculated using GetArrayManifold.m or can
% be calculated internally by passing an empty matrix [] as the argument
%
% SetupAudio is the Samples x Mics matrix of audio to be used for MVDR/LCMV
% weight calculation
%
% BFWeights and Evalstruct are struct arrays of dimensions methods x angles
% for a certain scene
% (any other variations should be taken care of in a higher level script)
%
% ArrayMan is returned as an optional third argument
%
% Rxx matrix may optionally be returned and passed into the function

if nargin<4
    getRxx=1;
else
    getRxx=0;
end

if isempty(ArrayMan)
    
    ArrayMan=GetArrayManifold(Config);
    
end

if getRxx
    if ismember('lcmv',Config.Filter.Method) || ismember('mvdr',Config.Filter.Method)
        Rxx=CalcRxx(SetupAudio,Config,ArrayMan.Setup);
    else
        Rxx=[];
    end
end

if ~iscell(Rxx)
    Rxx1=Rxx;
    Rxx=cell(size(Config.Scene.Targets.Angle,1));
    for itarget=1:size(Config.Scene.Targets.Angle,1)
        Rxx{itarget}=Rxx1;
    end
end

BFWeights = cell(length(Config.Filter.Method),size(Config.Scene.Targets.Angle,1));
EvalStruct = cell(length(Config.Filter.Method),size(Config.Scene.Targets.Angle,1));

Config.Filter.CurrInts = [];
for itarget=1:size(Config.Scene.Targets.Angle,1)
    Config.Filter.CurrTarget=Config.Scene.Targets.Angle(itarget,:);
    if ~isempty(Config.Scene.Interferers.Angle)
        Config.Filter.CurrInts=Config.Scene.Interferers.Angle(itarget,:);
    end

    ils=0;
    for imethod=1:length(Config.Filter.Method)
        bftype=Config.Filter.Method{imethod};
        if strcmp(bftype,'ls') || strcmp(bftype,'mn')
            ils = ils+1;
            Config.Filter.(sprintf(bftype)).CurrMode = Config.Filter.(sprintf(bftype)).mode{ils};
            Config.Filter.(sprintf(bftype)).CurrOrder = Config.Filter.(sprintf(bftype)).order(ils);
        end
        BFWeights{imethod,itarget}=BeamformerWeights(Rxx{itarget},Config,ArrayMan,bftype);
        EvalStruct{imethod,itarget}=BFWeightEvaluation(BFWeights{imethod,itarget},Config,ArrayMan.Monitor);
    end
end
end
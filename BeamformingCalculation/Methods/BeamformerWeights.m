function [BFWeights] = BeamformerWeights(Rxx,Config,ArrayMan,bftype)
% [BFAudio,BFWeights]=BeamformerWeights(Audio,Config,ArrayMan,bftype)
% Calculates beamformer weights and applies them to an audio signal.
%
% input arguments:  
%    Rxx: array covariance matrix of the input audio signal.
%    Config: configuration struct containing all the settings.
%    ArrayMan: struct containing the array manifold transfer functions.
%    bftype: string containing beamformer type to be applied. bftype 
%       options: 'ds', 'sda', 'mvdr', 'lcmv', 'acc', 'ls' and 'mn'.
%
% output arguments: 
%    BFWeights: beamforming weights struct. 
%
% Philip Coleman, CVSSP, 11/2014


switch bftype
    case 'ds'
        BFWeights=DelaySumNarrowband(Config,ArrayMan.Setup);
    case 'sda'
        BFWeights=SDANarrowband(Config,ArrayMan.Setup);
    case 'mvdr'
        Config.Filter.CurrInts=[];
        BFWeights=LCMVNarrowband(Config,ArrayMan.Setup,Rxx);
    case 'lcmv'
        BFWeights=LCMVNarrowband(Config,ArrayMan.Setup,Rxx);
    case 'acc'
        BFWeights=ACCNarrowband(Config,ArrayMan.Setup);
    case 'ls'
        BFWeights=LSNarrowband(Config,ArrayMan.Setup);
    case 'mn'
        BFWeights=MNNarrowband(Config,ArrayMan.Setup);
    %case 'maxsnr'
    %    BFWeights=MaxSNR(Config);
    otherwise
        error('invalid beamformer type specified');
end




end
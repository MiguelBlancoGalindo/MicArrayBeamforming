function [RespLogDiff,flog]=DirRespRefVsOtherRadius(BFEvalRef,BFEvalOtherR,ConfigOtherR,FreqAxis,RadiusSize)
% function [RespLogDiff,flog]=DirRespRefVsOtherRadius(BFEvalRef,BFEvalOtherR,ConfigOtherR,FreqAxis,RadiusSize)
% function that plots the difference between the response with an array of
% different radius than the reference case. 
%
% input arguments:
%   BFEvalRef: struct containing containing the directional response of the 
%       reference case.
%   BFEvalOtherR: struct containing containing the directional response of 
%       the other simulation set with different array radius.
%   ConfigOtherR: cell array containing all the configuration settings 
%       for the simulation set with different array radius.
%   FreqAxis: string taking values of 'linear' or 'log' depending on
%       whether BFEvalRef and BFEvalOtherR have a linear of logarithmic 
%       frequency axis. 
%   RadiusSize: string which can be set to 'half' or 'twice'. This is used 
%       to adjust the frequency axis of the reference response in order to 
%       shift it to match it to that of the array with different radius.
%
% output arguments:
%   RespLogDiff: directivity responses difference between reference and
%       simulation set with different array radius, after adjusting their
%       frequency ranges. 
%   flog: logarithmic frequency vector. 
%
% Note BFEvalRef and BFEvalOtherR are structs from the output of 
% rearrange_directional_response_measures.m. 

nangles = size(BFEvalRef,8);
f=ConfigOtherR{1,1,1}.Filter.FreqVector;
freqs = length(f);
foctaves = [32 63 125 250 500 1000 2000 4000 8000 16000];
fbandcutoff = [foctaves.*2^(-1/2) foctaves(end).*2^(1/2)];
%f=linspace(fbandcutoff(1),fbandcutoff(end),1000);  %test to check that
%with linear frequency vector starting and ending at fbandcutoff the number
if strcmp(FreqAxis,'linear')
    
    %of frequencies per octave band is constant.
    ifirstoctave = find(foctaves>=f(1),1,'first');
    ilastoctave = find(foctaves<=f(end),1,'last');
    foctaves = foctaves(ifirstoctave:ilastoctave);
    %iband = assign_freq_to_bands(f,fbandcutoff);
    fperoctave = 20;
    foct0 = fperoctave - find(f(1)<foctaves(1)/2*2.^(1/fperoctave:1/fperoctave:1),1,'first') + 1;
    foctinf = find(f(end)<foctaves(end)*2.^(1/fperoctave:1/fperoctave:1),1,'first')-1;
    flog = logspace(log10(f(1)),log10(f(end)),fperoctave*(length(foctaves)-1) + foct0 + foctinf);
    %flog = foctaves(floor(length(foctaves)/2)+1).*2.^(-floor(length(foctaves)/2)-(foct0-1)/fperoctave:1/fperoctave:length(foctaves)-floor(length(foctaves)/2)-1+foctinf/fperoctave);
    flog_prime = zeros(1,length(flog));
    iflog_prime = zeros(1,length(flog));
    for iflog=1:length(flog)
        [flog_prime(iflog),iflog_prime(iflog)] = min(abs(f-flog(iflog)));
    end

    %iflog_prime = unique(iflog_prime);
    RespRefLog = BFEvalRef.Respaz(:,:,:,:,:,:,:,iflog_prime);
    RespRefShifted = zeros(size(BFEvalRef.Respaz));
    RespOtherR = zeros(size(BFEvalOtherR.Respaz));
    iflog_HFcut = unique(iflog_prime(1:end-fperoctave));
    iflog_LFcut = unique(iflog_prime(fperoctave+1:end));
    clear iflog_prime

    if strcmp (RadiusSize, 'twice')
        RespRefLogShifted = BFEvalRef.Respaz(:,:,:,:,:,:,:,iflog_HFcut,:);
        RespOtherRLog = BFEvalOtherR.Respaz(:,:,:,:,:,:,:,iflog_LFcut,:);
        RespRefShifted(:,:,:,:,:,:,:,iflog_HFcut,:) = BFEvalRef.Respaz(:,:,:,:,:,:,:,iflog_HFcut,:);
        RespOtherR(:,:,:,:,:,:,:,iflog_LFcut,:) = BFEvalOtherR.Respaz(:,:,:,:,:,:,:,iflog_LFcut,:);
        %flog = flog(1:end-fperoctave);
        for iflog=1:length(iflog_HFcut)
            [flog_prime(iflog),iflog_prime(iflog)] = min(abs(f(iflog_HFcut(iflog))-flog));
        end
        flog = flog(iflog_prime);
    elseif strcmp (RadiusSize, 'half')
        RespRefLogShifted = BFEvalRef.Respaz(:,:,:,:,:,:,:,iflog_LFcut,:);
        RespOtherRLog = BFEvalOtherR.Respaz(:,:,:,:,:,:,:,iflog_HFcut,:);
        RespRefShifted(:,:,:,:,:,:,:,iflog_LFcut,:) = BFEvalRef.Respaz(:,:,:,:,:,:,:,iflog_LFcut,:);
        RespOtherR(:,:,:,:,:,:,:,iflog_HFcut,:) = BFEvalOtherR.Respaz(:,:,:,:,:,:,:,iflog_HFcut,:);
        %flog = flog(fperoctave+1:end);
        for iflog=1:length(iflog_LFcut)
            [flog_prime(iflog),iflog_prime(iflog)] = min(abs(f(iflog_LFcut(iflog))-flog));
        end
        flog = flog(iflog_prime);
    else
        error('choose between "twice" or "half" to set "RadiusSize"');
    end
    
    RespDiff = abs(RespRefShifted-RespOtherR);
    
elseif strcmp(FreqAxis,'log')
    fperoctave = length(f(f>=500 & f<1000));
    RespRefLog = BFEvalRef.Respaz;
    if strcmp (RadiusSize, 'twice')
        RespRefLogShifted = BFEvalRef.Respaz(:,:,:,:,:,:,:,fperoctave+1:end,:);
        RespOtherRLog = BFEvalOtherR.Respaz(:,:,:,:,:,:,:,1:end-fperoctave,:);
        flog = f(1:end-fperoctave);
    elseif strcmp (RadiusSize, 'half')
        RespRefLogShifted = BFEvalRef.Respaz(:,:,:,:,:,:,:,1:end-fperoctave,:);
        RespOtherRLog = BFEvalOtherR.Respaz(:,:,:,:,:,:,:,fperoctave+1:end,:);
        flog = f(fperoctave+1:end);
    else
        error('choose between "twice" or "half" to set "RadiusSize"');
    end
    
end

RespLogDiff = abs(RespRefLogShifted./RespOtherRLog);

end
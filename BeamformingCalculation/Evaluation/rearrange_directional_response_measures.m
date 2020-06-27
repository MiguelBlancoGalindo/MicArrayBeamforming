function BFEvalAllArrays2 = rearrange_directional_response_measures(ConfigAllArrays,BFEvalAllArrays)
% function BFEvalAllArrays2 = rearrange_directional_response_measures(ConfigAllArrays,BFEvalAllArrays)
% function that rearranges the evaluation cell array into a struct with
% each field being a metric with a multidimensional array. 
%
% input arguments: 
%   ConfigAllArrays: cell array containing all the configuration settings 
%       for each microphone array, beamforming method, etc.
%   BFEvalAllArrays: cell array containing several performance
%       metrics for each array, beamforming method, etc.
%
% output arguments:
%   BFEvalAllArrays2: struct containing the same information as 
%       BFEvalAllArrays but in a struct form with each metric being a 
%       multidimensional array. 

f=ConfigAllArrays{1,1,1,1}.Filter.FreqVector;
freqs = length(f); 
N=ndims(BFEvalAllArrays);
narray_geo = size(BFEvalAllArrays,1);
nM = size(ConfigAllArrays,2);
nr = size(ConfigAllArrays,3);
nepsilons = size(ConfigAllArrays,4);
nmethods = size(BFEvalAllArrays{1,1,1,1},1);
nsources = size(BFEvalAllArrays{1,1,1,1},2);
nnoises = length(BFEvalAllArrays{1,1,1,1}{1,1});
doaz=0; doel=0; do3d=0; doFR=0; doACRange=0;

if isfield(BFEvalAllArrays{1,1,1,1}{1,1}(1),'BWaz')
    nangles = size(BFEvalAllArrays{1,1,1,1}{1,1}(1).Respaz,2);
    DIaz = zeros(narray_geo,nM,nr,nmethods,nsources,nnoises,nepsilons,freqs); 
    BWaz = zeros(narray_geo,nM,nr,nmethods,nsources,nnoises,nepsilons,freqs); 
    BW3az = zeros(narray_geo,nM,nr,nmethods,nsources,nnoises,nepsilons,freqs); 
    BW15az = zeros(narray_geo,nM,nr,nmethods,nsources,nnoises,nepsilons,freqs); 
    SLaz = zeros(narray_geo,nM,nr,nmethods,nsources,nnoises,nepsilons,freqs); 
    WNG = zeros(narray_geo,nM,nr,nmethods,nsources,nnoises,nepsilons,freqs); 
    Contrastaz = zeros(narray_geo,nM,nr,nmethods,nsources,nnoises,nepsilons,freqs);
    Respaz = zeros(narray_geo,nM,nr,nmethods,nsources,nnoises,nepsilons,freqs,nangles);
    GLAz = zeros(narray_geo,nM,nr,nmethods,nsources,nnoises,nepsilons,freqs);
    GLAzf = zeros(narray_geo,nM,nr,nmethods,nsources,nnoises,nepsilons,freqs);
    if isfield(BFEvalAllArrays{1,1,1,1}{1,1}(1),'FRespaz')
        doFR=1;
        FRangeaz = zeros(narray_geo,nM,nr,nmethods,nsources,nnoises,nepsilons,2);
        FRespaz = zeros(narray_geo,nM,nr,nmethods,nsources,nnoises,nepsilons,freqs);
    end
    if isfield(BFEvalAllArrays{1,1,1,1}{1,1}(1),'ContrastRangeaz')
        ContrastRangeaz = zeros(narray_geo,nM,nr,nmethods,nsources,nnoises,nepsilons,freqs);
        doACRange=1;
    end
    doaz=1;
end
if isfield(BFEvalAllArrays{1,1,1,1}{1,1}(1),'BWel')
    nangles = size(BFEvalAllArrays{1,1,1,1}{1,1}(1).Respel,2);
    DIel = zeros(narray_geo,nM,nr,nmethods,nsources,nnoises,nepsilons,freqs); 
    BWel = zeros(narray_geo,nM,nr,nmethods,nsources,nnoises,nepsilons,freqs); 
    BW3el = zeros(narray_geo,nM,nr,nmethods,nsources,nnoises,nepsilons,freqs); 
    BW15el = zeros(narray_geo,nM,nr,nmethods,nsources,nnoises,nepsilons,freqs); 
    SLel = zeros(narray_geo,nM,nr,nmethods,nsources,nnoises,nepsilons,freqs); 
    Contrastel = zeros(narray_geo,nM,nr,nmethods,nsources,nnoises,nepsilons,freqs);
    WNG = zeros(narray_geo,nM,nr,nmethods,nsources,nnoises,nepsilons,freqs); 
    Respel = nan(narray_geo,nM,nr,nmethods,nsources,nnoises,nepsilons,freqs,nangles);
    GLEl = zeros(narray_geo,nM,nr,nmethods,nsources,nnoises,nepsilons,freqs);
    GLElf = zeros(narray_geo,nM,nr,nmethods,nsources,nnoises,nepsilons,freqs);
    if isfield(BFEvalAllArrays{1,1,1,1}{1,1}(1),'FRespel')
        doFR=1;
        FRangeel = zeros(narray_geo,nM,nr,nmethods,nsources,nnoises,nepsilons,2);
        FRespel = zeros(narray_geo,nM,nr,nmethods,nsources,nnoises,nepsilons,freqs); 
    end
    if isfield(BFEvalAllArrays{1,1,1,1}{1,1}(1),'ContrastRangeel')
        ContrastRangeel = zeros(narray_geo,nM,nr,nmethods,nsources,nnoises,nepsilons,freqs);
        doACRange=1;
    end
    doel=1;
end

if isfield(BFEvalAllArrays{1,1,1,1}{1,1}(1),'DI3d')
    nlook = size(BFEvalAllArrays{1,1,1,1}{1,1}(1).Resp3d,2);
    Resp3d = zeros(narray_geo,nM,nr,nmethods,nsources,nnoises,nepsilons,freqs,nlook);
    DI3d = zeros(narray_geo,nM,nr,nmethods,nsources,nnoises,nepsilons,freqs);
    WNG = zeros(narray_geo,nM,nr,nmethods,nsources,nnoises,nepsilons,freqs); 
    Contrast3d = zeros(narray_geo,nM,nr,nmethods,nsources,nnoises,nepsilons,freqs);
    do3d=1;
    if isfield(BFEvalAllArrays{1,1,1,1}{1,1}(1),'ContrastRange3d')
        ContrastRange3d = zeros(narray_geo,nM,nr,nmethods,nsources,nnoises,nepsilons,freqs);
        doACRange=1;
    end
end


for imethod=1:nmethods
    for isource=1:nsources
        for inoise=1:nnoises
            for iarray_geo=1:narray_geo
                for iM=1:nM
                    for ir=1:nr
                        for ie = 1:nepsilons
                            if doaz
                                DIaz(iarray_geo,iM,ir,imethod,isource,inoise,ie,:) = BFEvalAllArrays{iarray_geo,iM,ir,ie}{imethod,isource}(inoise).DIaz;
                                BWaz(iarray_geo,iM,ir,imethod,isource,inoise,ie,:) = BFEvalAllArrays{iarray_geo,iM,ir,ie}{imethod,isource}(inoise).BWaz;
                                BW3az(iarray_geo,iM,ir,imethod,isource,inoise,ie,:) = BFEvalAllArrays{iarray_geo,iM,ir,ie}{imethod,isource}(inoise).BW3az;
                                BW15az(iarray_geo,iM,ir,imethod,isource,inoise,ie,:) = BFEvalAllArrays{iarray_geo,iM,ir,ie}{imethod,isource}(inoise).BW15az;
                                SLaz(iarray_geo,iM,ir,imethod,isource,inoise,ie,:) = BFEvalAllArrays{iarray_geo,iM,ir,ie}{imethod,isource}(inoise).SLaz;
                                WNG(iarray_geo,iM,ir,imethod,isource,inoise,ie,:) = BFEvalAllArrays{iarray_geo,iM,ir,ie}{imethod,isource}(inoise).WNG;
                                Contrastaz(iarray_geo,iM,ir,imethod,isource,inoise,ie,:) = BFEvalAllArrays{iarray_geo,iM,ir,ie}{imethod,isource}(inoise).Contrastaz;
                                Respaz(iarray_geo,iM,ir,imethod,isource,inoise,ie,:,:) = BFEvalAllArrays{iarray_geo,iM,ir,ie}{imethod,isource}(inoise).Respaz;
                                if ~isempty(BFEvalAllArrays{iarray_geo,iM,ir,ie}{imethod,isource}(inoise).GL.Az)
                                    GLlength = length(BFEvalAllArrays{iarray_geo,iM,ir,ie}{imethod,isource}(inoise).GL.Az);
                                    GLAz(iarray_geo,iM,ir,imethod,isource,inoise,ie,1:GLlength) = BFEvalAllArrays{iarray_geo,iM,ir,ie}{imethod,isource}(inoise).GL.Az;
                                    GLAzf(iarray_geo,iM,ir,imethod,isource,inoise,ie,1:GLlength) = BFEvalAllArrays{iarray_geo,iM,ir,ie}{imethod,isource}(inoise).GL.Azf;
                                end
                                if doFR
                                    FRespaz(iarray_geo,iM,ir,imethod,isource,inoise,ie,:) = BFEvalAllArrays{iarray_geo,iM,ir,ie}{imethod,isource}(inoise).FRespaz;
                                    FRangeaz(iarray_geo,iM,ir,imethod,isource,inoise,ie,:) = BFEvalAllArrays{iarray_geo,iM,ir,ie}{imethod,isource}(inoise).FRangeaz; 
                                end
                                if doACRange; ContrastRangeaz(iarray_geo,iM,ir,imethod,isource,inoise,ie,:) = BFEvalAllArrays{iarray_geo,iM,ir,ie}{imethod,isource}(inoise).ContrastRangeaz; end
                            end
                            if doel
                                DIel(iarray_geo,iM,ir,imethod,isource,inoise,ie,:) = BFEvalAllArrays{iarray_geo,iM,ir,ie}{imethod,isource}(inoise).DIel;
                                BWel(iarray_geo,iM,ir,imethod,isource,inoise,ie,:) = BFEvalAllArrays{iarray_geo,iM,ir,ie}{imethod,isource}(inoise).BWel;
                                BW3el(iarray_geo,iM,ir,imethod,isource,inoise,ie,:) = BFEvalAllArrays{iarray_geo,iM,ir,ie}{imethod,isource}(inoise).BW3el;
                                BW15el(iarray_geo,iM,ir,imethod,isource,inoise,ie,:) = BFEvalAllArrays{iarray_geo,iM,ir,ie}{imethod,isource}(inoise).BW15el;
                                SLel(iarray_geo,iM,ir,imethod,isource,inoise,ie,:) = BFEvalAllArrays{iarray_geo,iM,ir,ie}{imethod,isource}(inoise).SLel;
                                WNG(iarray_geo,iM,ir,imethod,isource,inoise,ie,:) = BFEvalAllArrays{iarray_geo,iM,ir,ie}{imethod,isource}(inoise).WNG;
                                Contrastel(iarray_geo,iM,ir,imethod,isource,inoise,ie,:) = BFEvalAllArrays{iarray_geo,iM,ir,ie}{imethod,isource}(inoise).Contrastel;
                                Respel(iarray_geo,iM,ir,imethod,isource,inoise,ie,:,1:size(BFEvalAllArrays{iarray_geo,iM,ir,ie}{imethod,isource}(inoise).Respel,2)) = BFEvalAllArrays{iarray_geo,iM,ir,ie}{imethod,isource}(inoise).Respel;
                                if ~isempty(BFEvalAllArrays{iarray_geo,iM,ir,ie}{imethod,isource}(inoise).GL.El)
                                    GLlength = length(BFEvalAllArrays{iarray_geo,iM,ir,ie}{imethod,isource}(inoise).GL.El);
                                    GLEl(iarray_geo,iM,ir,imethod,isource,inoise,ie,1:GLlength) = BFEvalAllArrays{iarray_geo,iM,ir,ie}{imethod,isource}(inoise).GL.El;
                                    GLElf(iarray_geo,iM,ir,imethod,isource,inoise,ie,1:GLlength) = BFEvalAllArrays{iarray_geo,iM,ir,ie}{imethod,isource}(inoise).GL.Elf;
                                end
                                if doFR
                                    FRespel(iarray_geo,iM,ir,imethod,isource,inoise,ie,:) = BFEvalAllArrays{iarray_geo,iM,ir,ie}{imethod,isource}(inoise).FRespel;
                                    FRangeel(iarray_geo,iM,ir,imethod,isource,inoise,ie,:) = BFEvalAllArrays{iarray_geo,iM,ir,ie}{imethod,isource}(inoise).FRangeel; 
                                end
                                if doACRange; ContrastRangeel(iarray_geo,iM,ir,imethod,isource,inoise,ie,:) = BFEvalAllArrays{iarray_geo,iM,ir,ie}{imethod,isource}(inoise).ContrastRangeel; end
                            end
                            if do3d
                                DI3d(iarray_geo,iM,ir,imethod,isource,inoise,ie,:) = BFEvalAllArrays{iarray_geo,iM,ir,ie}{imethod,isource}(inoise).DI3d;
                                WNG(iarray_geo,iM,ir,imethod,isource,inoise,ie,:) = BFEvalAllArrays{iarray_geo,iM,ir,ie}{imethod,isource}(inoise).WNG; 
                                Contrast3d(iarray_geo,iM,ir,imethod,isource,inoise,ie,:) = BFEvalAllArrays{iarray_geo,iM,ir,ie}{imethod,isource}(inoise).Contrast3d;
                                if doACRange; ContrastRange3d(iarray_geo,iM,ir,imethod,isource,inoise,ie,:) = BFEvalAllArrays{iarray_geo,iM,ir,ie}{imethod,isource}(inoise).ContrastRange3d; end
                                Resp3d(iarray_geo,iM,ir,imethod,isource,inoise,ie,:,:) = BFEvalAllArrays{iarray_geo,iM,ir,ie}{imethod,isource}(inoise).Resp3d;
                            end

                        end
                    end
                end
            end
        end
    end
end

if doaz
    BFEvalAllArrays2.DIaz = DIaz;
    BFEvalAllArrays2.BWaz = BWaz;
    BFEvalAllArrays2.BW3az = BW3az;
    BFEvalAllArrays2.BW15az = BW15az;
    BFEvalAllArrays2.SLaz = SLaz;
    BFEvalAllArrays2.WNG = WNG;
    BFEvalAllArrays2.Contrastaz = Contrastaz;
    if doACRange; BFEvalAllArrays2.ContrastRangeaz = ContrastRangeaz; end
    BFEvalAllArrays2.Respaz = Respaz;
    BFEvalAllArrays2.GL.Az=GLAz;
    BFEvalAllArrays2.GL.Azf=GLAzf;
    if doFR
        BFEvalAllArrays2.FRespaz = FRespaz;
        BFEvalAllArrays2.FRangeaz = FRangeaz;
    end

elseif doel
    BFEvalAllArrays2.DIel = DIel;
    BFEvalAllArrays2.BWel = BWel;
    BFEvalAllArrays2.BW15el = BW15el;
    BFEvalAllArrays2.SLel = SLel;
    BFEvalAllArrays2.WNG = WNG;
    BFEvalAllArrays2.Contrastel = Contrastel;
    if doACRange; BFEvalAllArrays2.ContrastRangeel = ContrastRangeel; end
    BFEvalAllArrays2.Respel = Respel;
    BFEvalAllArrays2.GL.El=GLEl;
    BFEvalAllArrays2.GL.Elf=GLElf;
    if doFR
        BFEvalAllArrays2.FRespel = FRespel;
        BFEvalAllArrays2.FRangeel = FRangeel;
    end
end
if do3d
    BFEvalAllArrays2.Resp3d = Resp3d;
    BFEvalAllArrays2.DI3d = DI3d;
    BFEvalAllArrays2.WNG = WNG;
    BFEvalAllArrays2.Contrast3d = Contrast3d;
    if doACRange; BFEvalAllArrays2.ContrastRange3d = ContrastRange3d; end
end

end
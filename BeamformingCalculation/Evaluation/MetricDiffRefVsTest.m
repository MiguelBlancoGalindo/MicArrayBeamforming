function BFEvalAllArraysDiff = MetricDiffRefVsTest(BFEvalAllArraysTest, varargin)
%function that calculates the different between test response and reference
%response. The reference response can be included in BFEvalAllArraysTest as
%it is the case for non-ideal simulations where the first element of this
%array corresponds to reference case. Alternatively, and often used when 
%comparing array radii or number of microphones, the reference response can 
%be passed as an independent structure of 8D Arrays.

if isfield(BFEvalAllArraysTest,'DIaz')
    plane = 'az';
elseif isfield(BFEvalAllArraysTest,'DIel')
    plane = 'el';
end

%method to include part of a string as a field of a struct by using sprintf
%of the concatenated string ("DI"=part of variable + "plane"=newstring) and
%its output (string) include it in brackets to convert it to a field
narray_geo = size(BFEvalAllArraysTest.(sprintf(strcat('DI',plane))),1);
nM = size(BFEvalAllArraysTest.(sprintf(strcat('DI',plane))),2);
nr = size(BFEvalAllArraysTest.(sprintf(strcat('DI',plane))),3);
nmethods = size(BFEvalAllArraysTest.(sprintf(strcat('DI',plane))),4);
nsources = size(BFEvalAllArraysTest.(sprintf(strcat('DI',plane))),5);
nnoises = size(BFEvalAllArraysTest.(sprintf(strcat('DI',plane))),6);
nepsilons = size(BFEvalAllArraysTest.(sprintf(strcat('DI',plane))),7);
nfreqs = size(BFEvalAllArraysTest.(sprintf(strcat('DI',plane))),8);
if isfield(BFEvalAllArraysTest,'Respaz'); naz = size(BFEvalAllArraysTest.Respaz,9); end
if isfield(BFEvalAllArraysTest,'Respel'); nel = size(BFEvalAllArraysTest.Respaz,9); end

doaz=0;
doel=0;
do3d=0;
doFR=0;
if isfield(BFEvalAllArraysTest,'DIaz')
    doaz=1; 
    if isfield(BFEvalAllArraysTest,'FRespaz')
        doFR=1;
    end
end
if isfield(BFEvalAllArraysTest,'DIel')
    doel=1;
    if isfield(BFEvalAllArraysTest,'FRespel')
        doFR=1;
    end
end

if isfield(BFEvalAllArraysTest,'DIel')
    do3d=1;
end

% When passing a reference response struct, e.g. for array aperture or number of sensors 
if nargin==2
    BFEvalAllArraysRef = varargin{1};

if doaz && doel
    DI3d = BFEvalAllArraysTest.DI3d - BFEvalAllArraysRef.DI3d;
    WNG = BFEvalAllArraysTest.WNG - BFEvalAllArraysRef.WNG;
    Resp3d = BFEvalAllArraysTest.Resp3d - BFEvalAllArraysRef.Resp3d;
    Resp3D = BFEvalAllArraysTest.Resp3D - BFEvalAllArraysRef.Resp3D;
end
if doaz
    DIaz = BFEvalAllArraysTest.DIaz - BFEvalAllArraysRef.DIaz;
    BWaz = BFEvalAllArraysTest.BWaz - BFEvalAllArraysRef.BWaz;
    BW15az = BFEvalAllArraysTest.BW15az - BFEvalAllArraysRef.BW15az;
    SLaz = BFEvalAllArraysTest.SLaz - BFEvalAllArraysRef.SLaz;
    Contrastaz = BFEvalAllArraysTest.Contrastaz - BFEvalAllArraysRef.Contrastaz;
    Respaz = BFEvalAllArraysTest.Respaz - BFEvalAllArraysRef.Respaz;
    if doel==0;  WNG = BFEvalAllArraysTest.WNG - BFEvalAllArraysRef.WNG; end
    if doFR==1; FRespaz = BFEvalAllArraysTest.FRespaz - BFEvalAllArraysRef.FRespaz; end
end
if doel
    DIel = BFEvalAllArraysTest.DIel - BFEvalAllArraysRef.DIel;
    BWel = BFEvalAllArraysTest.BWel - BFEvalAllArraysRef.BWel;
    BW15el = BFEvalAllArraysTest.BW15el - BFEvalAllArraysRef.BW15el;
    SLel = BFEvalAllArraysTest.SLel - BFEvalAllArraysRef.SLel;
    Contrastel = BFEvalAllArraysTest.Contrastel - BFEvalAllArraysRef.Contrastel;
    Respel = BFEvalAllArraysTest.Respel - BFEvalAllArraysRef.Respel;
    if doaz==0;  WNG = BFEvalAllArraysTest.WNG - BFEvalAllArraysRef.WNG; end
    if doFR==1;  FRespel = BFEvalAllArraysTest.FRespel - BFEvalAllArraysRef.FRespel; end
end

% When using the same response struct, e.g. for non-ideal simulations 
elseif nargin==1 
    if doaz && doel
        DI3d = zeros(narray_geo,nM,nr,nmethods,nsources,nnoises-1,nepsilons,nfreqs);
        WNG = zeros(narray_geo,nM,nr,nmethods,nsources,nnoises-1,nepsilons,nfreqs);
        Resp3d = zeros(narray_geo,nM,nr,nmethods,nsources,nnoises-1,nepsilons,nfreqs,naz*nel);
        Resp3D = zeros(narray_geo,nM,nr,nmethods,nsources,nnoises-1,nepsilons,nfreqs,nel,naz);
        for inoise=1:nnoises-1
            DI3d(:,:,:,:,:,inoise,:,:) = BFEvalAllArraysTest.DI3d(:,:,:,:,:,inoise+1,:,:) - BFEvalAllArraysTest.DI3d(:,:,:,:,:,1,:,:);
            WNG(:,:,:,:,:,inoise,:,:) = BFEvalAllArraysTest.WNG(:,:,:,:,:,inoise+1,:,:) - BFEvalAllArraysTest.WNG(:,:,:,:,:,1,:,:);
            Resp3d(:,:,:,:,:,inoise,:,:) = BFEvalAllArraysTest.Resp3d(:,:,:,:,:,inoise+1,:,:) - BFEvalAllArraysTest.Resp3d(:,:,:,:,:,1,:,:);
            Resp3D(:,:,:,:,:,inoise,:,:) = BFEvalAllArraysTest.Resp3D(:,:,:,:,:,inoise+1,:,:) - BFEvalAllArraysTest.Resp3D(:,:,:,:,:,1,:,:);
        end
    end
    if doaz
        DIaz = zeros(narray_geo,nM,nr,nmethods,nsources,nnoises-1,nepsilons,nfreqs);
        BWaz = zeros(narray_geo,nM,nr,nmethods,nsources,nnoises-1,nepsilons,nfreqs);
        BW15az = zeros(narray_geo,nM,nr,nmethods,nsources,nnoises-1,nepsilons,nfreqs);
        SLaz = zeros(narray_geo,nM,nr,nmethods,nsources,nnoises-1,nepsilons,nfreqs);
        Contrastaz = zeros(narray_geo,nM,nr,nmethods,nsources,nnoises-1,nepsilons,nfreqs);
        Respaz = zeros(narray_geo,nM,nr,nmethods,nsources,nnoises-1,nepsilons,nfreqs,naz);
        if doel==0; WNG = zeros(narray_geo,nM,nr,nmethods,nsources,nnoises-1,nepsilons,nfreqs); end
        if doFR==1; FRespaz = zeros(narray_geo,nM,nr,nmethods,nsources,nnoises-1,nepsilons,nfreqs); end
        
        for inoise=1:nnoises-1
            DIaz(:,:,:,:,:,inoise,:,:) = BFEvalAllArraysTest.DIaz(:,:,:,:,:,inoise+1,:,:) - BFEvalAllArraysTest.DIaz(:,:,:,:,:,1,:,:);
            BWaz(:,:,:,:,:,inoise,:,:) = BFEvalAllArraysTest.BWaz(:,:,:,:,:,inoise+1,:,:) - BFEvalAllArraysTest.BWaz(:,:,:,:,:,1,:,:);
            BW15az(:,:,:,:,:,inoise,:,:) = BFEvalAllArraysTest.BW15az(:,:,:,:,:,inoise+1,:,:) - BFEvalAllArraysTest.BW15az(:,:,:,:,:,1,:,:);
            SLaz(:,:,:,:,:,inoise,:,:) = BFEvalAllArraysTest.SLaz(:,:,:,:,:,inoise+1,:,:) - BFEvalAllArraysTest.SLaz(:,:,:,:,:,1,:,:);
            Contrastaz(:,:,:,:,:,inoise,:,:) = BFEvalAllArraysTest.Contrastaz(:,:,:,:,:,inoise+1,:,:) - BFEvalAllArraysTest.Contrastaz(:,:,:,:,:,1,:,:);
            Respaz(:,:,:,:,:,inoise,:,:) = BFEvalAllArraysTest.Respaz(:,:,:,:,:,inoise+1,:,:) - BFEvalAllArraysTest.Respaz(:,:,:,:,:,1,:,:);
            if doel==0;  WNG(:,:,:,:,:,inoise,:,:) = BFEvalAllArraysTest.WNG(:,:,:,:,:,inoise+1,:,:) - BFEvalAllArraysTest.WNG(:,:,:,:,:,1,:,:); end
            if doFR==1; FRespaz(:,:,:,:,:,inoise,:,:) = BFEvalAllArraysTest.FRespaz(:,:,:,:,:,inoise+1,:,:) - BFEvalAllArraysTest.FRespaz(:,:,:,:,:,1,:,:); end
        end
    end
    if doel
        DIel = zeros(narray_geo,nM,nr,nmethods,nsources,nnoises-1,nepsilons,nfreqs);
        BWel = zeros(narray_geo,nM,nr,nmethods,nsources,nnoises-1,nepsilons,nfreqs);
        BW15el = zeros(narray_geo,nM,nr,nmethods,nsources,nnoises-1,nepsilons,nfreqs);
        SLel = zeros(narray_geo,nM,nr,nmethods,nsources,nnoises-1,nepsilons,nfreqs);
        Contrastel = zeros(narray_geo,nM,nr,nmethods,nsources,nnoises-1,nepsilons,nfreqs);
        Respel = zeros(narray_geo,nM,nr,nmethods,nsources,nnoises-1,nepsilons,nfreqs,nel);
        if doaz==0; WNG = zeros(narray_geo,nM,nr,nmethods,nsources,nnoises-1,nepsilons,nfreqs); end
        if doFR==1; FRespel = zeros(narray_geo,nM,nr,nmethods,nsources,nnoises-1,nepsilons,nfreqs); end
        
        for inoise=1:nnoises-1
            DIel(:,:,:,:,:,inoise,:,:) = BFEvalAllArraysTest.DIel(:,:,:,:,:,inoise+1,:,:) - BFEvalAllArraysTest.DIel(:,:,:,:,:,1,:,:);
            BWel(:,:,:,:,:,inoise,:,:) = BFEvalAllArraysTest.BWel(:,:,:,:,:,inoise+1,:,:) - BFEvalAllArraysTest.BWel(:,:,:,:,:,1,:,:);
            BW15el(:,:,:,:,:,inoise,:,:) = BFEvalAllArraysTest.BW15el(:,:,:,:,:,inoise+1,:,:) - BFEvalAllArraysTest.BW15el(:,:,:,:,:,1,:,:);
            SLel(:,:,:,:,:,inoise,:,:) = BFEvalAllArraysTest.SLel(:,:,:,:,:,inoise+1,:,:) - BFEvalAllArraysTest.SLel(:,:,:,:,:,1,:,:);
            Contrastel(:,:,:,:,:,inoise,:,:) = BFEvalAllArraysTest.Contrastel(:,:,:,:,:,inoise+1,:,:) - BFEvalAllArraysTest.Contrastel(:,:,:,:,:,1,:,:);
            Respel(:,:,:,:,:,inoise,:,:) = BFEvalAllArraysTest.Respel(:,:,:,:,:,inoise+1,:,:) - BFEvalAllArraysTest.Respel(:,:,:,:,:,1,:,:);
            if doaz==0;  WNG(:,:,:,:,:,inoise,:,:) = BFEvalAllArraysTest.WNG(:,:,:,:,:,inoise+1,:,:) - BFEvalAllArraysTestf.WNG(:,:,:,:,:,1,:,:); end
            if doFR==1;  FRespel(:,:,:,:,:,inoise,:,:) = BFEvalAllArraysTest.FRespel(:,:,:,:,:,inoise+1,:,:) - BFEvalAllArraysTest.FRespel(:,:,:,:,:,1,:,:); end
        end
    end
end


if doaz
    BFEvalAllArraysDiff.DIaz = DIaz;
    BFEvalAllArraysDiff.BWaz = BWaz;
    BFEvalAllArraysDiff.BW15az = BW15az;
    BFEvalAllArraysDiff.SLaz = SLaz;
    BFEvalAllArraysDiff.WNG = WNG;
    BFEvalAllArraysDiff.Contrastaz = Contrastaz;
    BFEvalAllArraysDiff.Respaz = Respaz;
    if doFR; BFEvalAllArraysDiff.FRespaz = FRespaz; end
end
if doel
    BFEvalAllArraysDiff.DIel = DIel;
    BFEvalAllArraysDiff.BWel = BWel;
    BFEvalAllArraysDiff.BW15el = BW15el;
    BFEvalAllArraysDiff.SLel = SLel;
    BFEvalAllArraysDiff.Contrastel = Contrastel;
    BFEvalAllArraysDiff.Respel = Respel;
    if doFR; BFEvalAllArraysDiff.FRespel = FRespel; end
end
if do3d
    BFEvalAllArraysDiff.DI3d = DI3d;
    BFEvalAllArraysDiff.Resp3d = Resp3d;
    BFEvalAllArraysDiff.Resp3D = Resp3D;
end

end



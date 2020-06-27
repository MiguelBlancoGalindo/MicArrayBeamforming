function y = BeamformerOutput(x,w)
% function y = BeamformerOutput(x,w)
% function that filters the input signal with the microphone array weights.
%
% input arguments:
%   x: microphone array signals.
%   w: beamforming weights.
%
% output argument:
%   y: beamformer output signal. 

y = sum(fftfilt(w,x),2);
y = 0.999.*y./max(abs(y));
end
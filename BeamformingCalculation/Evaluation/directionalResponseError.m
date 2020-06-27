function NSE = directionalResponseError(d,d_d)
% function NSE = directionalResponseError(d,d_d)
% function that calculates the normalised squared error in dB between the
% synthesised and the target directivity function.
%
% input arguments:
%   d: synthesised directivity response (with the beamformer weights).
%   d_d: target (or desired) directivity response.
%
% output arguments:
%   NSE: normalised squared error (in dB) between actual and target responses.


error=abs(d - abs(d_d));

inan = all(isnan(error),1);
NSE = 10*log10(sum(error(:,~inan).^2,2)./sum(abs(d(:,~inan).^2),2));

end
function WNG=getcurrWNG(w,a_target)

%WNG=10*log10(abs(w'*a_target)^2./(w'*w));
WNG=10*log10(mean(abs(w'*a_target).^2)./(w'*w));
function a = GetDiffMicCoeff(N,type)
%function that obtains the coefficients to synthesise a directivity pattern
%of type 'type' and order 'N' for spherical sound fields. They are of the form:
    %a(1) + a(2)*cos(alpha) + a(3)*cos.^2(alpha) + ... + a(N-1)*cos.^N(alpha)
%Coefficients obtained from "Fundamentals of Spherical Array Processing" by 
%Rafaely

if strcmp(type,'hyp')
    if N==0
        a=1;
    elseif N==1
        a=1/4.*[1 3];
    elseif N==2
        a=1/6.*[-1 2 5];
    elseif N==3
        a=1/32.*[-3 -15 15 35];
    elseif N==4
        a=1/40.*[3 -12 42 28 63];
    elseif N==5
        a=1/96.*[5 35 -70 -210 105 231];
    end
end
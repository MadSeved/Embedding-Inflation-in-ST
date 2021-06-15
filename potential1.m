% This function contains the KKLT potential used in "racetrack.m". The
% variables are the real and imaginary part of the K?hler modulus.

function V = potential1(X, Y)
global A
global a
global W0
global E
global alpha

V=E./(X.^alpha) + ((a.*A.*exp(-a.*X))./(2.*X.^2)).*(A.*exp(-a.*X).*((X.*a)./3+1)   + W0.*cos(a.*Y)) ;
end
% This function contains the Racetrack potential used in "racetrack.m". The
% variables are the real and imaginary part of the K?hler modulus.

function V = potential(X, Y)
global A
global a
global B
global b
global W0
global E
global alpha

V=E./(X.^alpha) + (exp(-a.*X))./(6.*(X.^2)).*(a.*(A.^2).*(a.*X+3).*exp(-a.*X)+3.*W0.*a.*A.*cos(a.*Y)) + (exp(-b.*X))./(6.*(X.^2)).*(b.*(B.^2).*(b.*X+3).*exp(-b.*X)+3.*W0.*b.*B.*cos(b.*Y))+  (exp(-(a+b).*X))./(6.*(X.^2)).*(A.*B.*(2.*a.*b.*X+3.*a+3.*b).*cos((a-b).*Y));
end
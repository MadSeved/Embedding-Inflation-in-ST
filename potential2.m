% This function contains the Racetrack potential used in "racetrack2.m". The
% variables are the real and imaginary part of the Kahler moduli.

function V = potential2(X1, X2, Y1, Y2)
global A
global a
global B
global b
global W0
global D

V = D./(((X2.^(3/2))-(X1.^(3/2))).^2) + (216./(((X2.^(3/2))-(X1.^(3/2))).^2)) .* ((B.^2).*b.*(b.*(X2.^2)+2.*b.*(X1.^(3/2)).*(X2.^(1/2))+3.*X2).*exp(-2.*b.*X2)      +     (A.^2).*a.*(a.*(X1.^2)+2.*a.*(X2.^(3/2)).*(X1.^(1/2))+3.*X1).*exp(-2.*a.*X1)         +        3.*B.*b.*W0.*X2.*exp(-b.*X2).*cos(b.*Y2)    +     3.*A.*a.*W0.*X1.*exp(-a.*X1).*cos(a.*Y1)      +      3.*A.*B.*exp(-a.*X1-b.*X2).*(a.*X1+b.*X2+2.*a.*b.*X1.*X2).*cos(-a.*Y1+b.*Y2));
end
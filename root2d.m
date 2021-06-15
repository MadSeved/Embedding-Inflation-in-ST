% This function contains the first order X derivatives of the Racetrack potential
% used in "racetrack2.m". The Y coordinates are fixed to y1 and y2 respectively.

function F = root2d(x)
global A
global a
global B
global b
global W0
global D
global y1
global y2

% Here we calculate the derivatives of the potential.
syms X1 X2
dX1V = diff(D./(((X2.^(3/2))-(X1.^(3/2))).^2) + (216./(((X2.^(3/2))-(X1.^(3/2))).^2)) .* ((B.^2).*b.*(b.*(X2.^2)+2.*b.*(X1.^(3/2)).*(X2.^(1/2))+3.*X2).*exp(-2.*b.*X2)      +     (A.^2).*a.*(a.*(X1.^2)+2.*a.*(X2.^(3/2)).*(X1.^(1/2))+3.*X1).*exp(-2.*a.*X1)         +        3.*B.*b.*W0.*X2.*exp(-b.*X2).*cos(b.*y2)    +     3.*A.*a.*W0.*X1.*exp(-a.*X1).*cos(a.*y1)      +      3.*A.*B.*exp(-a.*X1-b.*X2).*(a.*X1+b.*X2+2.*a.*b.*X1.*X2).*cos(-a.*y1+b.*y2)), X1);
dX2V = diff(D./(((X2.^(3/2))-(X1.^(3/2))).^2) + (216./(((X2.^(3/2))-(X1.^(3/2))).^2)) .* ((B.^2).*b.*(b.*(X2.^2)+2.*b.*(X1.^(3/2)).*(X2.^(1/2))+3.*X2).*exp(-2.*b.*X2)      +     (A.^2).*a.*(a.*(X1.^2)+2.*a.*(X2.^(3/2)).*(X1.^(1/2))+3.*X1).*exp(-2.*a.*X1)         +        3.*B.*b.*W0.*X2.*exp(-b.*X2).*cos(b.*y2)    +     3.*A.*a.*W0.*X1.*exp(-a.*X1).*cos(a.*y1)      +      3.*A.*B.*exp(-a.*X1-b.*X2).*(a.*X1+b.*X2+2.*a.*b.*X1.*X2).*cos(-a.*y1+b.*y2)), X2);

% Here we replace the symbolic variables with the input variables.
f(1)=subs(dX1V, [X1, X2], [x(1), x(2)]);
f(2)=subs(dX2V, [X1, X2], [x(1), x(2)]);

%Here we express the result above as a numerical value.
F(1)=double(f(1));
F(2)=double(f(2));
F=F*10^54; % We rescale "F" in order to make it easier for "fsolve".
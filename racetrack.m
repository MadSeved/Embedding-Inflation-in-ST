% This script plots the KKLT and Racetrack scalar potentials against the
% inflaton, i.e. the K?hler modulus, for a certain set of parameters. Then
% it calculates the slow-roll parameters for the Racetrack potential at its
% saddle point. The scripts also calculates the mass matrix for the
% Racetrack potential.

% Line 8 and 9 removes variables from the Workspace and closes all Matlab windows.
clear
close all

% Line 12 to 26 defines the set of parameters and makes them global.
global A
global a
global B
global b
global W0
global E
global alpha

A=1/50; 
a=(2*pi)/100;
B=-35/1000;
b=(2*pi)/90;
W0=-1/25000;
E=4.14668*10^(-12);
alpha=2;

[X, Y]=meshgrid(100:1:160, -30:1:30); % This line defines the meshgrid that we plot the Racetrack potential against. X is the real part of the K?hler potential, Y is the imaginary part.
V=potential(X, Y)*10^16; % We define "V" to be the Racetrack potential. The factor of 10^16 is there for scaling purposes.
figure
surf(X, Y, V) % In this line we make the actual surface plot for the Racetrack potential.
title('Racetrack potential')
zlabel('V')
xlabel('X')
ylabel('Y')

V1=potential1(X, Y)*10^14; % We define "V1" to be the KKLT potential. The factor of 10^14 is there for scaling purposes.
figure
surf(X, Y, V1) % In this line we make the surface plot for the KKLT potential.
title('KKLT potential')
zlabel('V')
xlabel('X')
ylabel('Y')

% In line 48 to 54 we define the first and second order derivatives of the
% Racetrack potential. The variables are the real and imaginary parts of
% the K?hler potential.
syms X Y
VderX=diff((E./(X.^alpha) + (exp(-a.*X))./(6.*(X.^2)).*(a.*(A.^2).*(a.*X+3).*exp(-a.*X)+3.*W0.*a.*A.*cos(a.*Y)) + (exp(-b.*X))./(6.*(X.^2)).*(b.*(B.^2).*(b.*X+3).*exp(-b.*X)+3.*W0.*b.*B.*cos(b.*Y))+  (exp(-(a+b).*X))./(6.*(X.^2)).*(A.*B.*(2.*a.*b.*X+3.*a+3.*b).*cos((a-b).*Y))), X);
VderY=diff((E./(X.^alpha) + (exp(-a.*X))./(6.*(X.^2)).*(a.*(A.^2).*(a.*X+3).*exp(-a.*X)+3.*W0.*a.*A.*cos(a.*Y)) + (exp(-b.*X))./(6.*(X.^2)).*(b.*(B.^2).*(b.*X+3).*exp(-b.*X)+3.*W0.*b.*B.*cos(b.*Y))+  (exp(-(a+b).*X))./(6.*(X.^2)).*(A.*B.*(2.*a.*b.*X+3.*a+3.*b).*cos((a-b).*Y))), Y);
Vderder11X=diff(VderX, X);
Vderder11Y=diff(VderX, Y);
Vderder12X=diff(VderY, X);
Vderder12Y=diff(VderY, Y);

% In line 57 and 58 we set X and Y to be the position of the saddle point.
X=123.22;
Y=0;

% "Vderder" is defined as the mass matrix of the Racetrack potential. In
% line 62 to 65 we define its elements at the saddle point.
Vderder(1,1)=subs(Vderder11X);
Vderder(1,2)=subs(Vderder11Y);
Vderder(2,1)=subs(Vderder12X);
Vderder(2,2)=subs(Vderder12Y);

% In line 68 and 69 we express the mass matrix and its eigenvalues.
Q=double(Vderder)
eig(Q)

% In line 73 and 74 we express the value of the Racetrack potential at the
% saddle point.
derp=['The value of the racetrack potential at the saddlepoint is ', num2str(potential(X, Y))];
disp(derp)

% In line 78 to 81 we calculate the slow-roll parameters at the saddle point and
% display their values.
epsilonRT=((X^2)/3)*(subs(VderY)/potential(X, Y))^2;
etaRT=((2*(X^2))/3)*(Vderder(2,2)/potential(X, Y));
epsilon_saddle=vpa(epsilonRT)
eta_saddle=vpa(etaRT)
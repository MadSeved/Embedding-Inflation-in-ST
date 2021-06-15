% This script finds the X-plane coordinates of the minimum and saddle
% point of the Racetrack potential with two Kahler moduli
% as inflaton fields for a certain set of parameters and fixed Y-plane
% coordinates. X and Y are the real and imaginary part of the inflaton respectively.
% Then the script calculate the slow-roll parameters
% for a potential with four variables (the variables being the real and
% imaginary parts of the Kahler moduli). Then the script plot the
% potential against the imaginary part of the inflaton whilst fixing the
% real part to its minimum, plot the potential against the real part of the
% inflaton whilst fixing the imaginary part to one of its minumum and also
% plot the potential against the real part of the inflaton whilst fixing
% the imaginart part to one of its maximum.

% Line 15 and 16 removes variables from the Workspace and closes all Matlab windows.
clear
close all

% Here we define the set of parameters and make them global. We also
% introduce y1 and y2 as global constants. They will play the roles of
% fixed Y coordinates when necessary.
global A
global a
global B
global b
global W0
global D
global y1
global y2
A=0.56;
a=(2*pi)/40;
B=7.46666*10^-5;
b=(2*pi)/258;
W0=5.22666*10^-6;
D=6.21019*10^-9;

% Here we fix the Y coordinates to a minimum.
y1=0;
y2=pi/b;

% Here we calculate and display the minimum in the X-plane for Y fixed at
% the minimum. This corresponds to a four dimensional minimum.
fun = @root2d; % fun contains the derivatives of the potential. 
x0 = [100,180]; % x0 acts as the starting guess for the Matlab solver "fsolve".
x = fsolve(fun,x0); % This matlab solver finds the minimum. It does so by finding where the derivatives are zero.
resultat=['The four dimensional minimum is placed at X1=', num2str(x(1)), ', X2=', num2str(x(2)),', Y1=0, Y2=129.'];
disp(resultat);

% Here we fix the Y coordinates to a maximum.
y1=pi/a;
y2=pi/b;

% Here we calculate and display the minimum in the X-plane for Y fixed at
% the maximum. This corresponds to a four dimensional saddle point.
fun = @root2d; % fun contains the derivatives of the potential.
x0 = [100,180]; % x0 acts as the starting guess for the Matlab solver "fsolve".
x = fsolve(fun,x0); % This matlab solver finds the minimum. It does so by finding where the derivatives are zero.
resultat=['The four dimensional saddle point is placed at X1=', num2str(x(1)), ', X2=', num2str(x(2)),', Y1=20, Y2=129.'];
disp(resultat);

% In line 61 and 62 we define the Kahler potential.
syms T1 T2 barT1 barT2
K=log(1296)-2*log((T2+barT2)^(3/2)   -    (T1+barT1)^(3/2));

% Here we define the elements of the Kahler metric.
barT1K=diff(K, barT1);
barT2K=diff(K, barT2);
K11=(diff(barT1K, T1));
K12=(diff(barT2K, T1));
K21=(diff(barT1K, T2));
K22=(diff(barT2K, T2));

% Here we construct the Kahler metric and its inverse "QQinv".
QQ(1,1)=K11;
QQ(1,2)=K12;
QQ(2,1)=K21;
QQ(2,2)=K22;
QQinv=inv(QQ);

% Here we introduce the Racetrack potential and then we change its
% variables from X, Y to T, barT.
syms X1 X2 Y1 Y2 T1 T2 barT1 barT2
V=D./(((X2.^(3/2))-(X1.^(3/2))).^2) + (216./(((X2.^(3/2))-(X1.^(3/2))).^2)) .* ((B.^2).*b.*(b.*(X2.^2)+2.*b.*(X1.^(3/2)).*(X2.^(1/2))+3.*X2).*exp(-2.*b.*X2)      +     (A.^2).*a.*(a.*(X1.^2)+2.*a.*(X2.^(3/2)).*(X1.^(1/2))+3.*X1).*exp(-2.*a.*X1)         +        3.*B.*b.*W0.*X2.*exp(-b.*X2).*cos(b.*Y2)    +     3.*A.*a.*W0.*X1.*exp(-a.*X1).*cos(a.*Y1)      +      3.*A.*B.*exp(-a.*X1-b.*X2).*(a.*X1+b.*X2+2.*a.*b.*X1.*X2).*cos(-a.*Y1+b.*Y2));
nyV=subs(V, [X1, X2, Y1, Y2], [(T1+barT1)/2, (T2+barT2)/2, (T1-barT1)/(2*sqrt(-1)), (T2-barT2)/(2*sqrt(-1))]);

% Here we calculate epsilon.
epsilon=((QQinv(1,1)*diff(nyV, T1)*diff(nyV, barT1)   +   QQinv(1,2)*diff(nyV, T1)*diff(nyV, barT2)    +   QQinv(2,1)*diff(nyV, T2)*diff(nyV, barT1)    +     QQinv(2,2)*diff(nyV, T2)*diff(nyV, barT2))/nyV^2);

% In line 90 to 108 we calculate the elements in the matrix from which
% \eta is obtained.
N11=(QQinv(1,1)*diff(diff(nyV, T1), barT1)    +    QQinv(1,2)*diff(diff(nyV, T1), barT2))/(nyV);
N12=(QQinv(1,1)*diff(diff(nyV, T2), barT1)    +    QQinv(1,2)*diff(diff(nyV, T2), barT2))/(nyV);
N21=(QQinv(2,1)*diff(diff(nyV, T1), barT1)    +    QQinv(2,2)*diff(diff(nyV, T1), barT2))/(nyV);
N22=(QQinv(2,1)*diff(diff(nyV, T2), barT1)    +    QQinv(2,2)*diff(diff(nyV, T2), barT2))/(nyV);

barN11=(QQinv(1,1)*diff(diff(nyV, barT1), T1)    +    QQinv(2,1)*diff(diff(nyV, barT1), T2))/(nyV);
barN12=(QQinv(1,1)*diff(diff(nyV, barT2), T1)    +    QQinv(2,1)*diff(diff(nyV, barT2), T2))/(nyV);
barN21=(QQinv(1,2)*diff(diff(nyV, barT1), T1)    +    QQinv(2,2)*diff(diff(nyV, barT1), T2))/(nyV);
barN22=(QQinv(1,2)*diff(diff(nyV, barT2), T1)    +    QQinv(2,2)*diff(diff(nyV, barT2), T2))/(nyV);

N1bar1=(QQinv(1,1)*diff(diff(nyV, barT1), barT1) + QQinv(1,2)*diff(diff(nyV, barT2), barT1)   -   QQinv(1,1)*QQinv(1,1)*diff(diff(diff(K, T1), barT1), barT1)*diff(nyV, barT1)   -   QQinv(1,2)*QQinv(1,1)*diff(diff(diff(K, T1), barT2), barT1)*diff(nyV, barT1)    -   QQinv(1,1)*QQinv(2,1)*diff(diff(diff(K, T2), barT1), barT1)*diff(nyV, barT1)    -   QQinv(1,2)*QQinv(2,1)*diff(diff(diff(K, T2), barT2), barT1)*diff(nyV, barT1)    -   QQinv(1,1)*QQinv(1,2)*diff(diff(diff(K, T1), barT1), barT1)*diff(nyV, barT2)   -   QQinv(1,2)*QQinv(1,2)*diff(diff(diff(K, T1), barT2), barT1)*diff(nyV, barT2)    -   QQinv(1,1)*QQinv(2,2)*diff(diff(diff(K, T2), barT1), barT1)*diff(nyV, barT2)    -   QQinv(1,2)*QQinv(2,2)*diff(diff(diff(K, T2), barT2), barT1)*diff(nyV, barT2))/nyV;
N1bar2=(QQinv(1,1)*diff(diff(nyV, barT1), barT2) + QQinv(1,2)*diff(diff(nyV, barT2), barT2)   -   QQinv(1,1)*QQinv(1,1)*diff(diff(diff(K, T1), barT1), barT2)*diff(nyV, barT1)   -   QQinv(1,2)*QQinv(1,1)*diff(diff(diff(K, T1), barT2), barT2)*diff(nyV, barT1)    -   QQinv(1,1)*QQinv(2,1)*diff(diff(diff(K, T2), barT1), barT2)*diff(nyV, barT1)    -   QQinv(1,2)*QQinv(2,1)*diff(diff(diff(K, T2), barT2), barT2)*diff(nyV, barT1)    -   QQinv(1,1)*QQinv(1,2)*diff(diff(diff(K, T1), barT1), barT2)*diff(nyV, barT2)   -   QQinv(1,2)*QQinv(1,2)*diff(diff(diff(K, T1), barT2), barT2)*diff(nyV, barT2)    -   QQinv(1,1)*QQinv(2,2)*diff(diff(diff(K, T2), barT1), barT2)*diff(nyV, barT2)    -   QQinv(1,2)*QQinv(2,2)*diff(diff(diff(K, T2), barT2), barT2)*diff(nyV, barT2))/nyV;
N2bar1=(QQinv(2,1)*diff(diff(nyV, barT1), barT1) + QQinv(2,2)*diff(diff(nyV, barT2), barT1)   -   QQinv(2,1)*QQinv(1,1)*diff(diff(diff(K, T1), barT1), barT1)*diff(nyV, barT1)   -   QQinv(2,2)*QQinv(1,1)*diff(diff(diff(K, T1), barT2), barT1)*diff(nyV, barT1)    -   QQinv(2,1)*QQinv(2,1)*diff(diff(diff(K, T2), barT1), barT1)*diff(nyV, barT1)    -   QQinv(2,2)*QQinv(2,1)*diff(diff(diff(K, T2), barT2), barT1)*diff(nyV, barT1)    -   QQinv(2,1)*QQinv(1,2)*diff(diff(diff(K, T1), barT1), barT1)*diff(nyV, barT2)   -   QQinv(2,2)*QQinv(1,2)*diff(diff(diff(K, T1), barT2), barT1)*diff(nyV, barT2)    -   QQinv(2,1)*QQinv(2,2)*diff(diff(diff(K, T2), barT1), barT1)*diff(nyV, barT2)    -   QQinv(2,2)*QQinv(2,2)*diff(diff(diff(K, T2), barT2), barT1)*diff(nyV, barT2))/nyV;
N2bar2=(QQinv(2,1)*diff(diff(nyV, barT1), barT2) + QQinv(2,2)*diff(diff(nyV, barT2), barT2)   -   QQinv(2,1)*QQinv(1,1)*diff(diff(diff(K, T1), barT1), barT2)*diff(nyV, barT1)   -   QQinv(2,2)*QQinv(1,1)*diff(diff(diff(K, T1), barT2), barT2)*diff(nyV, barT1)    -   QQinv(2,1)*QQinv(2,1)*diff(diff(diff(K, T2), barT1), barT2)*diff(nyV, barT1)    -   QQinv(2,2)*QQinv(2,1)*diff(diff(diff(K, T2), barT2), barT2)*diff(nyV, barT1)    -   QQinv(2,1)*QQinv(1,2)*diff(diff(diff(K, T1), barT1), barT2)*diff(nyV, barT2)   -   QQinv(2,2)*QQinv(1,2)*diff(diff(diff(K, T1), barT2), barT2)*diff(nyV, barT2)    -   QQinv(2,1)*QQinv(2,2)*diff(diff(diff(K, T2), barT1), barT2)*diff(nyV, barT2)    -   QQinv(2,2)*QQinv(2,2)*diff(diff(diff(K, T2), barT2), barT2)*diff(nyV, barT2))/nyV;

barN1bar1=(QQinv(1,1)*diff(diff(nyV, T1), T1) + QQinv(2,1)*diff(diff(nyV, T2), T1)   -   QQinv(1,1)*QQinv(1,1)*diff(diff(diff(K, barT1), T1), T1)*diff(nyV, T1)   -   QQinv(2,1)*QQinv(1,1)*diff(diff(diff(K, barT1), T2), T1)*diff(nyV, T1)    -   QQinv(1,1)*QQinv(1,2)*diff(diff(diff(K, barT2), T1), T1)*diff(nyV, T1)    -   QQinv(2,1)*QQinv(1,2)*diff(diff(diff(K, barT2), T2), T1)*diff(nyV, T1)    -   QQinv(1,1)*QQinv(2,1)*diff(diff(diff(K, barT1), T1), T1)*diff(nyV, T2)   -   QQinv(2,1)*QQinv(2,1)*diff(diff(diff(K, barT1), T2), T1)*diff(nyV, T2)    -   QQinv(1,1)*QQinv(2,2)*diff(diff(diff(K, barT2), T1), T1)*diff(nyV, T2)    -   QQinv(2,1)*QQinv(2,2)*diff(diff(diff(K, barT2), T2), T1)*diff(nyV, T2))/nyV;
barN1bar2=(QQinv(1,1)*diff(diff(nyV, T1), T2) + QQinv(2,1)*diff(diff(nyV, T2), T2)   -   QQinv(1,1)*QQinv(1,1)*diff(diff(diff(K, barT1), T1), T2)*diff(nyV, T1)   -   QQinv(2,1)*QQinv(1,1)*diff(diff(diff(K, barT1), T2), T2)*diff(nyV, T1)    -   QQinv(1,1)*QQinv(1,2)*diff(diff(diff(K, barT2), T1), T2)*diff(nyV, T1)    -   QQinv(2,1)*QQinv(1,2)*diff(diff(diff(K, barT2), T2), T2)*diff(nyV, T1)    -   QQinv(1,1)*QQinv(2,1)*diff(diff(diff(K, barT1), T1), T2)*diff(nyV, T2)   -   QQinv(2,1)*QQinv(2,1)*diff(diff(diff(K, barT1), T2), T2)*diff(nyV, T2)    -   QQinv(1,1)*QQinv(2,2)*diff(diff(diff(K, barT2), T1), T2)*diff(nyV, T2)    -   QQinv(2,1)*QQinv(2,2)*diff(diff(diff(K, barT2), T2), T2)*diff(nyV, T2))/nyV;
barN2bar1=(QQinv(1,2)*diff(diff(nyV, T1), T1) + QQinv(2,2)*diff(diff(nyV, T2), T1)   -   QQinv(1,2)*QQinv(1,1)*diff(diff(diff(K, barT1), T1), T1)*diff(nyV, T1)   -   QQinv(2,2)*QQinv(1,1)*diff(diff(diff(K, barT1), T2), T1)*diff(nyV, T1)    -   QQinv(1,2)*QQinv(1,2)*diff(diff(diff(K, barT2), T1), T1)*diff(nyV, T1)    -   QQinv(2,2)*QQinv(1,2)*diff(diff(diff(K, barT2), T2), T1)*diff(nyV, T1)    -   QQinv(1,2)*QQinv(2,1)*diff(diff(diff(K, barT1), T1), T1)*diff(nyV, T2)   -   QQinv(2,2)*QQinv(2,1)*diff(diff(diff(K, barT1), T2), T1)*diff(nyV, T2)    -   QQinv(1,2)*QQinv(2,2)*diff(diff(diff(K, barT2), T1), T1)*diff(nyV, T2)    -   QQinv(2,2)*QQinv(2,2)*diff(diff(diff(K, barT2), T2), T1)*diff(nyV, T2))/nyV;
barN2bar2=(QQinv(1,2)*diff(diff(nyV, T1), T2) + QQinv(2,2)*diff(diff(nyV, T2), T2)   -   QQinv(1,2)*QQinv(1,1)*diff(diff(diff(K, barT1), T1), T2)*diff(nyV, T1)   -   QQinv(2,2)*QQinv(1,1)*diff(diff(diff(K, barT1), T2), T2)*diff(nyV, T1)    -   QQinv(1,2)*QQinv(1,2)*diff(diff(diff(K, barT2), T1), T2)*diff(nyV, T1)    -   QQinv(2,2)*QQinv(1,2)*diff(diff(diff(K, barT2), T2), T2)*diff(nyV, T1)    -   QQinv(1,2)*QQinv(2,1)*diff(diff(diff(K, barT1), T1), T2)*diff(nyV, T2)   -   QQinv(2,2)*QQinv(2,1)*diff(diff(diff(K, barT1), T2), T2)*diff(nyV, T2)    -   QQinv(1,2)*QQinv(2,2)*diff(diff(diff(K, barT2), T1), T2)*diff(nyV, T2)    -   QQinv(2,2)*QQinv(2,2)*diff(diff(diff(K, barT2), T2), T2)*diff(nyV, T2))/nyV;

% Here we change variables from T, barT to X, Y.
epsilon=subs(epsilon, [T1, T2, barT1, barT2], [X1+sqrt(-1)*Y1, X2+sqrt(-1)*Y2, X1-sqrt(-1)*Y1, X2-sqrt(-1)*Y2]);
N11=subs(N11, [T1, T2, barT1, barT2], [X1+sqrt(-1)*Y1, X2+sqrt(-1)*Y2, X1-sqrt(-1)*Y1, X2-sqrt(-1)*Y2]);
N12=subs(N12,[T1, T2, barT1, barT2], [X1+sqrt(-1)*Y1, X2+sqrt(-1)*Y2, X1-sqrt(-1)*Y1, X2-sqrt(-1)*Y2]);
N21=subs(N21,[T1, T2, barT1, barT2], [X1+sqrt(-1)*Y1, X2+sqrt(-1)*Y2, X1-sqrt(-1)*Y1, X2-sqrt(-1)*Y2]);
N22=subs(N22,[T1, T2, barT1, barT2], [X1+sqrt(-1)*Y1, X2+sqrt(-1)*Y2, X1-sqrt(-1)*Y1, X2-sqrt(-1)*Y2]);
barN11=subs(barN11,[T1, T2, barT1, barT2], [X1+sqrt(-1)*Y1, X2+sqrt(-1)*Y2, X1-sqrt(-1)*Y1, X2-sqrt(-1)*Y2]);
barN12=subs(barN12,[T1, T2, barT1, barT2], [X1+sqrt(-1)*Y1, X2+sqrt(-1)*Y2, X1-sqrt(-1)*Y1, X2-sqrt(-1)*Y2]);
barN21=subs(barN21,[T1, T2, barT1, barT2], [X1+sqrt(-1)*Y1, X2+sqrt(-1)*Y2, X1-sqrt(-1)*Y1, X2-sqrt(-1)*Y2]);
barN22=subs(barN22,[T1, T2, barT1, barT2], [X1+sqrt(-1)*Y1, X2+sqrt(-1)*Y2, X1-sqrt(-1)*Y1, X2-sqrt(-1)*Y2]);
N1bar1=subs(N1bar1,[T1, T2, barT1, barT2], [X1+sqrt(-1)*Y1, X2+sqrt(-1)*Y2, X1-sqrt(-1)*Y1, X2-sqrt(-1)*Y2]);
N1bar2=subs(N1bar2,[T1, T2, barT1, barT2], [X1+sqrt(-1)*Y1, X2+sqrt(-1)*Y2, X1-sqrt(-1)*Y1, X2-sqrt(-1)*Y2]);
N2bar1=subs(N2bar1,[T1, T2, barT1, barT2], [X1+sqrt(-1)*Y1, X2+sqrt(-1)*Y2, X1-sqrt(-1)*Y1, X2-sqrt(-1)*Y2]);
N2bar2=subs(N2bar2,[T1, T2, barT1, barT2], [X1+sqrt(-1)*Y1, X2+sqrt(-1)*Y2, X1-sqrt(-1)*Y1, X2-sqrt(-1)*Y2]);
barN1bar1=subs(barN1bar1,[T1, T2, barT1, barT2], [X1+sqrt(-1)*Y1, X2+sqrt(-1)*Y2, X1-sqrt(-1)*Y1, X2-sqrt(-1)*Y2]);
barN1bar2=subs(barN1bar2,[T1, T2, barT1, barT2], [X1+sqrt(-1)*Y1, X2+sqrt(-1)*Y2, X1-sqrt(-1)*Y1, X2-sqrt(-1)*Y2]);
barN2bar1=subs(barN2bar1,[T1, T2, barT1, barT2], [X1+sqrt(-1)*Y1, X2+sqrt(-1)*Y2, X1-sqrt(-1)*Y1, X2-sqrt(-1)*Y2]);
barN2bar2=subs(barN2bar2,[T1, T2, barT1, barT2], [X1+sqrt(-1)*Y1, X2+sqrt(-1)*Y2, X1-sqrt(-1)*Y1, X2-sqrt(-1)*Y2]);

% Here we prepare the matrix from which we obtain \eta,Q, as well as fix the
% coordinates X1, X2, Y1, Y2 to correspond to the four dimensional saddle
% point.
Q=zeros(4);
X1=x(1);
X2=x(2);
Y1=y1;
Y2=y2;

% Here we calculate and display the value of the slow-roll parameter
% \epsilon at the saddle point.
epsilon=subs(epsilon);
epsilon=double(epsilon)

% Here we evaluate the matrix Q at the saddle point as well as calculate
% the slow-roll parameter \eta as the same point.
Q(1,1)=subs(N11);
Q(1,2)=subs(N12);
Q(2,1)=subs(N21);
Q(2,2)=subs(N22);
Q(3,3)=subs(barN11);
Q(3,4)=subs(barN12);
Q(4,3)=subs(barN21);
Q(4,4)=subs(barN22);
Q(1,3)=subs(N1bar1);
Q(1,4)=subs(N1bar2);
Q(2,3)=subs(N2bar1);
Q(2,4)=subs(N2bar2);
Q(3,1)=subs(barN1bar1); 
Q(3,2)=subs(barN1bar2);
Q(4,1)=subs(barN2bar1);
Q(4,2)=subs(barN2bar2);
eig(Q);
eta=min(eig(Q))

% Here we display the value of the potential at the saddle point.
derp=['The value of the racetrack potential at the saddle point is ', num2str(potential2(X1, X2, Y1, Y2))];
disp(derp)

% Here we define the meshgrid for which we plot V2 against. V2 is the potential with Y fixed to a minimum.
[X1, X2]=meshgrid(90:1:170, 150:1:350);
V2=potential2(X1, X2, 0, 129)*10^14; % The factor of 10^14 is there for scaling purposes.

% Here we define the meshgrid for which we plot V1 against. V1 is the potential with X fixed to its minimum.
[Y1, Y2]=meshgrid(-40:2:40, 0:2:300);
V1=potential2(98.75839, 171.06117, Y1, Y2)*10^14; % The factor of 10^14 is there for scaling purposes.

% Line 177 to 184 alters some elements in V2 in order to get a surface plot
% with better colour scaling.
[M, N]=size(V2);
for i=1:1:M
    for j=1:1:N
        if V2(i,j)>0.51
V2(i,j)=0.51;
        end
    end
end
figure
surf(Y1, Y2, V1) % Here we make the V1 plot.
title('Plot of scalar potential versus Im(T1), Im(T2)')
zlabel('V')
xlabel('Y1')
ylabel('Y2')

figure
surf(X1, X2, V2) % Here we make the V2 plot.
axis([90 170 150 350 0 0.5]) % We compress the V-axis.
title('Plot of scalar potential versus Re(T1), Re(T2)')
zlabel('V')
xlabel('X1')
ylabel('X2')

figure
surf(X1, X2, V2) % Here we make a zoomed in version of the previous plot.
axis([90 130 150 250 0 0.05])
title('Plot of scalar potential versus Re(T1), Re(T2) zoomed in at the minimum')
zlabel('V')
xlabel('X1')
ylabel('X2')

% Here we define V2 as the potential with Y fixed to a maximum.
V2=potential2(X1, X2, 20, 129)*10^14;

% Line 214 to 221 alters some elements in V2 in order to get a surface plot
% with better colour scaling.
[M, N]=size(V2);
for i=1:1:M
    for j=1:1:N
        if V2(i,j)>0.51
V2(i,j)=0.51;
        end
    end
end
figure
surf(X1, X2, V2) % Here we make the new V2 plot.
axis([90 130 150 250 0 0.5])
title('Plot of scalar potential versus Re(T1), Re(T2)')
zlabel('V')
xlabel('X1')
ylabel('X2')

figure
surf(X1, X2, V2) % Here we make a zoomed in version of the previous plot.
axis([100 120 150 250 0.03 0.05])
title('Plot of scalar potential versus Re(T1), Re(T2) zoomed in at the minimum')
zlabel('V')
xlabel('X1')
ylabel('X2')
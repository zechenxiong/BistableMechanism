%% 1.test convergence or divergence
syms s A1 A2 b1 % s=0 correspond to the free end
A1 =0;
A2 =1;
b1 =4;
% phi = sqrt(s)*( A1*besselj(1/4,b1/2*s^2) + A2*besselj(-1/4,b1/2*s^2) );
phi = sqrt(s)*( A1*besselj(1/4,b1/2*s^2) + A2*besselj(-1/4,b1/2*s^2) );

% phi = besselj(1/4,s);
y = diff(phi,s);

vpa(subs(y,s,0.00001))
% our boundary condition that both end with zero rotation yields A2=0
% original encastre b.c. yields that A1=0
% and phi = sqrt(s)* A1*besselj(+-1/4,b1/2*s^2)

%% 2.find the zero of both
syms x
fplot(besselj(1/4,x),[0 10])
hold on

 myfun = @(x) besselj(1/4,x);
fzero(myfun,2); %ans = 2.7809;

%% 3.calculate Pcr from two conditions
clear
theta = 10; %deg

l = 20+7.5*sin(deg2rad(theta)); %(or 32.5?)
l2 = 75; % 
% l3 = 82.5;
L = 75+12.5;%-h/4; %(or 75 or 82.5) now its total ribbon length 87.5
E = 2000; % mm, N, MPa, s
v = 0.34; 
h = 15;
b = 0.381;

% method 1: from the lateral torsional buckling of the base short beam,
% which is pure bending
Pcr1 = pi/2/ (l2*cos(deg2rad(theta))) /l * sqrt( E * b^3*h/12 * E/2/(1+v) * b^3*h/3 );

% method 2: from the lateral torsional buckling of the limb long beam,
% first zero being 2.7809
Pcr2 = 5.5618/L^2 * sqrt( E * b^3*h/12 * E/2/(1+v) * b^3*h/3 );% / cos(deg2rad(theta));
RF = Pcr2/ cos(deg2rad(theta));% assuming the mode is still lateral buckling
% method 2: encastred beam lateral buckling
Pcr22 = 4.0126/L^2 * sqrt( E * b^3*h/12 * E/2/(1+v) * b^3*h/3 );% / cos(deg2rad(theta));

% method 3: column compressive buckling, i.e., when limb angle =90 deg
Pcr3 = pi^2/(1*L)^2 * E * b^3*h/12 / sin(deg2rad(theta));

%% 4.2021Apr16 plot UR1 theoretic along ribbon
% our boundary condition that both end with zero rotation yields A2=0 
% when s = L, i.e., s= L-z = L, z = 0, b1/2*s^2=2.7809
% calculating for phi(z), phi(s)

h=15;
Le = 75+12.5;%75+12.5-h/4; % equivalent ribbon length
L = 75+12.5;%-h/4; %(or 75 or 82.5) now its total ribbon length 87.5

E = 2000; % mm, N, MPa, s
v = 0.34; 
h = 15;
b = 0.381; %using Pcr2 with equivalent length after riveting
Pcr2 = 5.5618/Le^2 * sqrt( E * b^3*h/12 * E/2/(1+v) * b^3*h/3 );% Wrong!Using L instead of Le

b1 = 2.78088772399498*2/Le^2;
A1 = 0.100543107969855; %where to calculate A1, from 7?
theta = 10; %deg
syms z
s = Le-z;
Pcr2 = 5.5618/L^2 * sqrt( E * b^3*h/12 * E/2/(1+v) * b^3*h/3 );% / cos(deg2rad(theta));

phi =@(z) sqrt(Le-z).*A1.*besselj(1/4,b1/2.*(Le-z)^2);%.*cos(deg2rad(theta));
h = fplot(phi,[0 Le]);
x1 = h.XData;
y1 = h.YData;

%% 5.vpasolve integral, wrong, cannot solve A1
clear
% example:
% f1 = @(x_l) tan(asin(alpha.^2.*(a_l./4.*x_l-0.5.*x_l.^2)));
% f2 = @(x_l) integral(f1,0,x_l);
% [x1,y1] = fplot(f2,[0 a_l/2]);
% plot(x1,y1)
Le = 83.750000000000000;
b1 = 7.929498774782802e-04;
theta = 10;%deg
E = 2000; % mm, N, MPa, s
h = 15;%mm
v = 0.34; 

b = 0.381; %mm Pcr2 calculated from above-above section 0.1227

H = 12.5+75*sin(deg2rad(theta));
vB = H/cos(deg2rad(theta));
Pcr2 = 5.5618/Le^2 * sqrt( E * b^3*h/12 * E/2/(1+v) * b^3*h/3 );

syms A1 f1 s

f1 =  int( s*cos(sqrt(s)*A1* besselj(1/4,b1/2*s^2) ) ,s,Le,s);

eq1 = int(f1,s,Le,0) == -vB* E*h^3*b/12 /Pcr2;

A1 = double(vpasolve(eq1, A1, 0.1));

phi =@(z) sqrt(Le-z).*A1.*besselj(1/4,b1/2.*(Le-z)^2);
j = fplot(phi,[0 Le]);
x1 = j.XData;
y1 = j.YData;
% or [x1,y1]=fplot(phi,[0 Le])
%% 6.inverse calculation
clear
Le = 83.750000000000000;
b1 = 7.929498774782802e-04;
theta = 10;%deg
E = 2000; % mm, N, MPa, s
h = 15;%mm
b = 0.381; %mm Pcr2 calculated from above-above section 0.1227

H = 12.5+75*sin(deg2rad(theta));
vB = H/cos(deg2rad(theta));

i=1;
f3 = zeros(1,11);
for A1 = 0:0.1:1

    syms s a% A1
    f1 = int( s*cos(sqrt(s)*A1* besselj(1/4,b1/2*s^2) ) ,s,Le,a);
    % f1 = int( s*A1 ,s,Le,a);

    f2 = int(f1,a,Le,0);

    f3(i) = double(f2);
    i=i+1;
end
%     f3 =@(x) vpa(subs(f2,A1,x));%vpa() not working
% f3 = @(A1) f2;
% fplot(f3)


%% 7.calculating the A1 from energy method
clear
h = 15;%mm
Le = 87.5;%83.750000000000000; %beam length should be 12.5+75=87.5
E = 2000; % mm, N, MPa, s
b = 0.381; %mm Pcr2 calculated from above-above section 0.1227
v = 0.34; 
theta = 10;%deg

L = 75+12.5;%-h/4; %(or 75 or 82.5) now its total ribbon length 87.5, a way to increase A1
Pcr2 = 2.78088772399498*2/L^2 * sqrt( E * b^3*h/12 * E/2/(1+v) * b^3*h/3 );
b1 = 2.78088772399498*2/Le^2;


H = 12.5+75*sin(deg2rad(theta));
vB = H;%./cos(deg2rad(theta));  %a way to increase A1

syms A1 z
phi = sqrt(Le-z).*A1.*besselj(1/4,b1/2.*(Le-z)^2);
U = Pcr2^2/2/E/ (b^3*h/12) *int( sin(phi).^2.*(Le-z).^2 ,z,0,Le) + E/2/(1+v)*b^3*h/3 /2*int( diff(sin(phi),z,1).^2 ,z,0,Le);
%U = Pcr2^2/2/E/ (b^3*h/12) *int( phi.^2.*(Le-z).^2 ,z,0,Le) + E/2/(1+v)*b^3*h/3 /2*int( diff(phi,z,1).^2 ,z,0,Le);
% U = Pcr2^2/2/E/ (b^3*h/12) *int( phi.^2.*(Le-z).^2 ,z,0,Le) + E/2/(1+v)*b^3*h/3 /2*int( diff(phi,z,1).^2 ,z,0,Le) + Pcr2.^2/2/E/(h^3*b/12)*int((Le-z)^2,z,0,Le);
%U = Pcr2^2/2/E/ (b^3*h/12) *int( phi.^2.*(Le-z).^2 ,z,0,Le) + E/2/(1+v)*b^3*h/3 /2*int( diff(phi,z,1).^2 ,z,0,Le) + Pcr2.^2/2/E/(h^3*b/12)*int((Le-z)^2,z,0,Le);
%U = Pcr2^2/2/E/ (b^3*h/12) *int( phi.^2.*(L-z).^2 ,z,0,L) + E/2/(1+v)*b^3*h/3 /2*int( diff(phi,z,1).^2 ,z,0,L) + Pcr2.^2/2/E/(h^3*b/12)*int((L-z)^2,z,0,L);

T = Pcr2*vB;

eq1 = U == T;
A1 = double(vpasolve(eq1, A1, 0.100332869113339))
%% 8.calculating rotation? bending angle from A1 = 0.091582547512036~0.959~0.995~0.100332869113339~0.100543107969855
clear
h = 15;%mm
Le = 87.5;%83.750000000000000;%to change du_dz profile
E = 2000; % mm, N, MPa, s
b = 0.381; %mm Pcr2 calculated from above-above section 0.1227
v = 0.34; 
theta = 10;%deg

L = 75+12.5;%-h/4; %(or 75 or 82.5) now its total ribbon length 87.5, only influence Pcr2
Pcr2 = 5.5618/L^2 * sqrt( E * b^3*h/12 * E/2/(1+v) * b^3*h/3 );

b1 = 2.7809*2/Le^2;
A1 = 0.095288634708020;%0.0952886112145570;%0.096074615334912;%0.100543107969855;

syms z %a
phi = sqrt(Le-z).*A1.*besselj(1/4,b1/2.*(Le-z)^2);
du_dz = @(a) Pcr2/E/(b^3*h/12) *int(phi*(Le-z),z, 0, a);% * cos(subs(phi,z,a)); what is the point of adding cos(phi)

k = fplot(du_dz,[0 Le]);
x2(:,1) = k.XData;
x2(:,2) = k.YData;
%% 9.when a = 5 why there is negative value, Le should be 83.75, for that b1 is calculated from Le=83.75
clear
h = 15;%mm
Le = 83.750000000000000;%to change du_dz profile
E = 2000; % mm, N, MPa, s
b = 0.381; %mm Pcr2 calculated from above-above section 0.1227
v = 0.34; 
theta = 10;%deg

L = 75+12.5-h/4; %(or 75 or 82.5) now its total ribbon length 87.5, only influence Pcr2
Pcr2 = 5.5618/L^2 * sqrt( E * b^3*h/12 * E/2/(1+v) * b^3*h/3 );
b1 = 2.7809*2/Le^2;
A1 = 0.091582547512036;

a = 5;
syms z
phi = sqrt(Le-z).*A1.*besselj(1/4,b1/2.*(Le-z)^2);
du_dz = Pcr2/E/(b^3*h/12) *int(phi*(Le-z),z, 0, a);
%% 10.calculating the max energy state for snap-through buckling
% The bifurcation point of the snap-through buckling occurs when energy
% reaches the maximum. We can assume that the maximum of energy is achieved
% approximately when phi(z=20)= 0. The snap-through of the central ribbon lead to the
% subsequent snap-through of th limbs. Which has to be proved by fast
% filming camera
clear
h = 15;%mm
Le = 83.750000000000000;
E = 2000; % mm, N, MPa, s
b = 0.381; %mm Pcr2 calculated from above-above section 0.1227
v = 0.34; 
theta = 10;%deg

L = 75+12.5-h/4; %(or 75 or 82.5) now its total ribbon length 87.5, only influence Pcr2

Le_m = 12.5+75-h/4-20;
Pcr_m = 5.5618/Le_m^2 * sqrt( E * b^3*h/12 * E/2/(1+v) * b^3*h/3 );
Pcr2 = 5.5618/L^2 * sqrt( E * b^3*h/12 * E/2/(1+v) * b^3*h/3 );

b1_m = 2.7809*2/Le_m^2;


H = 12.5+75*sin(deg2rad(theta));
vB_m = H./cos(deg2rad(theta));  %a way to increase A1
vB = vB_m;

syms A3 z
phi_m = sqrt(Le_m-z).*A3.*besselj(1/4,b1_m/2.*(Le_m-z)^2);
U_m = Pcr_m^2/2/E/ (b^3*h/12) *int( sin(phi_m).^2.*(Le_m-z).^2 ,z,0,Le_m)...
      + E/2/(1+v)*b^3*h/3 /2*int( diff(sin(phi_m),z,1).^2 ,z,0,Le_m);
T_m = Pcr_m*vB_m;

eq3 = U_m == T_m;
A3 = double(vpasolve(eq3, A3, 0.1337));

syms z a
phi_m = sqrt(Le_m-z).*A3.*besselj(1/4,b1_m/2.*(Le_m-z)^2);
du_dz_m = @(a) Pcr_m/E/(b^3*h/12) *int(sin(phi_m)*(Le_m-z),z, 0, a) ...
        * cos(subs(phi_m,z,a));

p = fplot(du_dz_m,[0 Le_m]);
x3(:,1) = p.XData;
x3(:,2) = p.YData;
%% 11.iteration between du_dz and phi, following 8
% syms A1 z
H = 12.5+75*sin(deg2rad(theta));
vB_m = H./cos(deg2rad(theta));  %a way to increase A1
vB = vB_m;

% phi = sqrt(Le-z).*A1.*besselj(1/4,b1/2.*(Le-z)^2);
% U = Pcr2^2/2/E/ (b^3*h/12) *int( sin(phi).^2.*(Le-z).^2 ,z,0,Le) ...
%     + E/2/(1+v)*b^3*h/3 /2*int( diff(sin(phi),z,1).^2 ,z,0,Le);
% T = Pcr2*vB;
% eq1 = U == T;
% A1 = double(vpasolve(eq1, A1, 0.100332869113339));

syms z a
phi = @(z) sqrt(Le-z).*A1.*besselj(1/4,b1/2.*(Le-z)^2)...
      * cos(du_dz(Le)-du_dz(z));
j = fplot(phi,[0 Le]);
x4(:,1) = j.XData;
x4(:,2) = j.YData;

%% 12.calculating energy barrier, energy stored, and ratio contour
clear
theta = 0:3:45;
g_s = 2:0.5:10;

n = 1;
for i  = 1:length(theta)
    for j = 1:length(g_s)
%         [U_s_,D_U_,D_U_2,U_m_,D_U_3] = shellBarrEner(theta(i),g_s(j));
        [U_s_,D_U_] = shellBarrEner(theta(i),g_s(j));
        
        v(n,1) = theta(i);
        v(n,2) = g_s(j);
        v(n,3) = U_s_;
        v(n,4) = D_U_;
        
        n = n + 1;
    end
end
w = v((v(:,4)>=0),:);%choose the rows with non-negative values
%% 13.Calculating internal energy to compare to FEM
h = 15;
Le = 83.750000000000000;
E = 2000; % mm, N, MPa, s
b = 0.381; %mm Pcr2 calculated from above-above section 0.1227
v = 0.34; 
theta = 10;%deg

L = 75+12.5-h/4; %(or 75 or 82.5) now its total ribbon length 87.5, only influence Pcr2

Le_m = 12.5+75-h/4-20;

[U_s_,D_U_] = shellBarrEner(10,75/12.5);
U_s = U_s_*E* b^3*h/12 /12.5;

%% 14. following 8, recalculate bending angle, still A1 = 0.091582547512036~0.995~0.100332869113339~0.100543107969855, wrong, A1 not recalculated
clear
h = 15;%mm
Le = 87.5;%to change du_dz profile
E = 1738; % mm, N, MPa, s
b = 0.381; %mm Pcr2 calculated from above-above section 0.1227
v = 0.43; %this can also be changed
theta = 10;%deg
% A1 = 0.100543107969855;%changable
A1 = 0.095288634708020;%0.0952886112145570;%0.096074615334912;%0.100543107969855;


L = 75+12.5; %(or 75 or 82.5) now its total ribbon length 87.5, only influence Pcr2
Pcr2 = 5.5618/L^2 * sqrt( E * b^3*h/12 * E/2/(1+v) * b^3*h/3 );

b1 = 2.7809*2/Le^2;


syms z a
phi = sqrt(Le-z).*A1.*besselj(1/4,b1/2.*(Le-z)^2);
du_dz = @(a) Pcr2/E/(b^3*h/12) *int(phi*(Le-z),z, 0, a) * cos(subs(phi,z,a));

k = fplot(du_dz,[0 Le]);
x14(:,1) = k.XData;
x14(:,2) = k.YData;
%% 15. the whole process of calculating phi and bending angle psi, from 7+4+8
%7
clear

gamma = 6;
L1 = 12.5;
Le = L1*(gamma+1);%83.750000000000000; %beam length should be 12.5+75=87.5
v = 0.34;%0.43; 

h = 15*10;%mm
E = 1738*10;%1738; % mm, N, MPa, s
b = 0.381*10; %mm Pcr2 calculated from above-above section 0.1227

theta0 = [-9 -6 -4.5 -3 0 5 10 15 20 25 30 35 40 45];%deg
du_dz0 = zeros(length(theta0),2);


for i = 1: length(theta0)
    theta=theta0(i);
    L = Le;
    Pcr2 = 2.78088772399498*2/L^2 * sqrt( E * b^3*h/12 * E/2/(1+v) * b^3*h/3 );
    b1 = 2.78088772399498*2/Le^2;

    H = L1+L1*gamma*sin(deg2rad(theta));

    theta2= asin(1/gamma);
    vB = L1*gamma*(deg2rad(theta)+asin(1/gamma));  %a way to increase A1

    syms A1 z
    phi = sqrt(Le-z).*A1.*besselj(1/4,b1/2.*(Le-z)^2);
    U = Pcr2^2/2/E/ (b^3*h/12) *int( sin(phi).^2.*(Le-z).^2 ,z,0,Le) + E/2/(1+v)*b^3*h/3 /2*int( diff(sin(phi),z,1).^2 ,z,0,Le);%sin(phi) to increase psi?

    T = Pcr2*vB;

    eq1 = U == T;
    A1 = double(vpasolve(eq1, A1, 0.100332869113339))  % print A1 value
    du_dz0(i,2)=A1;
    %4
    Pcr2 = 5.5618/L^2 * sqrt( E * b^3*h/12 * E/2/(1+v) * b^3*h/3 );%

    phi =@(z) sqrt(Le-z).*A1.*besselj(1/4,b1/2.*(Le-z)^2);%

    % f = fplot(phi,[0 Le]);
    % x1 = f.XData;
    % y1 = f.YData;
    %8
    syms z
    phi1 = sqrt(Le-z).*A1.*besselj(1/4,b1/2.*(Le-z)^2);
    du_dz = @(a) Pcr2/E/(b^3*h/12) *int(phi1*(Le-z),z, 0, a);% * cos(subs(phi,z,a)); what is the point of adding cos(phi)

    du_dz0(i,1) = double(du_dz(Le));
end

%% 16. copied from 15, change back to phi', diff(phi,z,1)
clear

gamma = 6;

L1 = 12.5;%no influence on psi
Le = L1*(gamma+1);%83.750000000000000; %beam length should be 12.5+75=87.5
v = 0.3;%0.43; 

h = 15*10;%mm
E = 1738*10;%1738; % mm, N, MPa, s
b = 0.381*10; %mm Pcr2 calculated from above-above section 0.1227

theta0 = [-9 -6 -4.5 -3 0 5 10 15 20 25 30 35 40 45];%deg
du_dz0 = zeros(length(theta0),2);


for i = 1: length(theta0)
    theta=theta0(i);
    L = Le;
    Pcr2 = 2.78088772399498*2/L^2 * sqrt( E * b^3*h/12 * E/2/(1+v) * b^3*h/3 );
    b1 = 2.78088772399498*2/Le^2;

%     H = 12.5+75*sin(deg2rad(theta)); 
    H = L1+L1*gamma*sin(deg2rad(theta));

    theta2= asin(1/gamma);
    vB = L1*gamma*(deg2rad(theta)+asin(1/gamma));  %a way to increase A1

    syms A1 z
    phi = sqrt(Le-z).*A1.*besselj(1/4,b1/2.*(Le-z)^2);
    U = Pcr2^2/2/E/ (b^3*h/12) *int( sin(phi).^2.*(Le-z).^2 ,z,0,Le) + E/2/(1+v)*b^3*h/3 /2*int( diff(phi,z,1).^2 ,z,0,Le);%sin(phi) to increase psi?

    T = Pcr2*vB;

    eq1 = U == T;
    A1 = double(vpasolve(eq1, A1, 0.100332869113339))  % print A1 value
    du_dz0(i,2)=A1;
    %4
    Pcr2 = 5.5618/L^2 * sqrt( E * b^3*h/12 * E/2/(1+v) * b^3*h/3 );%

    phi =@(z) sqrt(Le-z).*A1.*besselj(1/4,b1/2.*(Le-z)^2);%

    % f = fplot(phi,[0 Le]);
    % x1 = f.XData;
    % y1 = f.YData;
    %8
    syms z
    phi1 = sqrt(Le-z).*A1.*besselj(1/4,b1/2.*(Le-z)^2);
    du_dz = @(a) Pcr2/E/(b^3*h/12) *int(phi1*(Le-z),z, 0, a);% * cos(subs(phi,z,a)); what is the point of adding cos(phi)

    du_dz0(i,1) = double(du_dz(Le));
end

%% 17. copied from 16, calculate psi with gamma
clear

gamma0 = [1.5 2 4 6 10];
L1 = 12.5;
v = 0.3;%0.43; 

h = 15;%mm
E = 1738;%1738; % mm, N, MPa, s
b = 0.381; %mm Pcr2 calculated from above-above section 0.1227

theta0 = 10;%deg
du_dz0 = zeros(length(gamma0),2);


for i = 1: length(gamma0)
    theta=theta0;
    gamma = gamma0(i);
    
    Le = L1*(gamma+1);
    L = Le;
    
    Pcr2 = 2.78088772399498*2/L^2 * sqrt( E * b^3*h/12 * E/2/(1+v) * b^3*h/3 );
    b1 = 2.78088772399498*2/Le^2;

    theta2= asin(1/gamma);
    vB = L1*gamma*(deg2rad(theta)+asin(1/gamma));  %a way to increase A1

    syms A1 z
    phi = sqrt(Le-z).*A1.*besselj(1/4,b1/2.*(Le-z)^2);
    U = Pcr2^2/2/E/ (b^3*h/12) *int( sin(phi).^2.*(Le-z).^2 ,z,0,Le) + E/2/(1+v)*b^3*h/3 /2*int( diff(phi,z,1).^2 ,z,0,Le);%sin(phi) to increase psi?

    T = Pcr2*vB;

    eq1 = U == T;
    A1 = double(vpasolve(eq1, A1, 0.100332869113339))  % print A1 value
    du_dz0(i,2)=A1;
    %4
    Pcr2 = 5.5618/L^2 * sqrt( E * b^3*h/12 * E/2/(1+v) * b^3*h/3 );%

    phi =@(z) sqrt(Le-z).*A1.*besselj(1/4,b1/2.*(Le-z)^2);%

    % f = fplot(phi,[0 Le]);
    % x1 = f.XData;
    % y1 = f.YData;
    %8
    syms z
    phi1 = sqrt(Le-z).*A1.*besselj(1/4,b1/2.*(Le-z)^2);
    du_dz = @(a) Pcr2/E/(b^3*h/12) *int(phi1*(Le-z),z, 0, a);% * cos(subs(phi,z,a)); what is the point of adding cos(phi)

    du_dz0(i,1) = double(du_dz(Le));
end

%% 18. copied from 16 plotting coutouring surface of psi w.r.t gamma and phi
clear

theta0 = [-4.5 -3 0 5 10 15 20 25 30 35 40 45];%deg
gamma0 = [1.5 2 4 6 8 10];

L1 = 12.5;%no influence on psi
v = 0.34;%0.43; 

h = 15*10;%mm
E = 1738*10;%1738; % mm, N, MPa, s
b = 0.381*10; %mm Pcr2 calculated from above-above section 0.1227


du_dz0 = zeros(length(theta0),length(gamma0));


for i = 1: length(theta0)
    for j = 1:length(gamma0)
        
        theta=theta0(i);
        gamma = gamma0(j);
        Le = L1*(gamma+1);%83.750000000000000; %beam length should be 12.5+75=87.5

        L = Le;
        Pcr2 = 2.78088772399498*2/L^2 * sqrt( E * b^3*h/12 * E/2/(1+v) * b^3*h/3 );
        b1 = 2.78088772399498*2/Le^2;

        theta2= asin(1/gamma);
        vB = L1*gamma*(deg2rad(theta)+asin(1/gamma))  %a way to increase A1

        syms A1 z
        phi = sqrt(Le-z).*A1.*besselj(1/4,b1/2.*(Le-z)^2);
        U = Pcr2^2/2/E/ (b^3*h/12) *int( sin(phi).^2.*(Le-z).^2 ,z,0,Le) + E/2/(1+v)*b^3*h/3 /2*int( diff(phi,z,1).^2 ,z,0,Le);%sin(phi) to increase psi?

        T = Pcr2*vB;

        eq1 = U == T;
        A1 = double(vpasolve(eq1, A1, 0.100332869113339))  % print A1 value
%         du_dz0(i,2)=A1;
        %4
        Pcr2 = 5.5618/L^2 * sqrt( E * b^3*h/12 * E/2/(1+v) * b^3*h/3 );%

        phi =@(z) sqrt(Le-z).*A1.*besselj(1/4,b1/2.*(Le-z)^2);%

        % f = fplot(phi,[0 Le]);
        % x1 = f.XData;
        % y1 = f.YData;
        %8
        syms z
        phi1 = sqrt(Le-z).*A1.*besselj(1/4,b1/2.*(Le-z)^2);
        du_dz = @(a) Pcr2/E/(b^3*h/12) *int(phi1*(Le-z),z, 0, a);% * cos(subs(phi,z,a)); what is the point of adding cos(phi)

        du_dz0(i,j) = double(du_dz(Le));
    end
end
%% 19. copied from 18 plotting coutouring surface of Ebarr w.r.t gamma and phi
clear

theta0 = [-4.5 -3 0 5 10 15 20 25 30 35 40 45];%deg
gamma0 = [1.5 2 4 6 8 10];

L1 = 12.5;%no influence on psi
v = 0.34;%0.43; 

h = 15*10;%mm
E = 2000*10;%1738; % mm, N, MPa, s
b = 0.381*10; %mm Pcr2 calculated from above-above section 0.1227


Ebarr = zeros(length(theta0),length(gamma0));


for i = 1: length(theta0)
    for j = 1:length(gamma0)
        
        theta=theta0(i);
        gamma = gamma0(j);
        Le = L1*(gamma+1);%83.750000000000000; %beam length should be 12.5+75=87.5
        L = Le;
        
        Ebarr(i,j)=shellBarrEner(theta,gamma);
        
    end
end

%% 20. prove psi and E_barr to be irrelavant to h/L1 and t/L1, copied from 18
clear

theta0 = [2 ];%deg
gamma0 = [5]; 

L1 = 12.5;%no influence on psi
v = 0.34;%0.43; 

h = 15*10;%mm
E = 1738*10;%1738; % mm, N, MPa, s
b = 0.381*2; %mm Pcr2 calculated from above-above section 0.1227


du_dz0 = zeros(length(theta0),length(gamma0));
Ebarr0 = zeros(length(theta0),length(gamma0));

for i = 1: length(theta0)
    for j = 1:length(gamma0)
        
        theta=theta0(i);
        gamma = gamma0(j);
        Le = L1*(gamma+1);%83.750000000000000; %beam length should be 12.5+75=87.5

        L = Le;
        Pcr2 = 2.78088772399498*2/L^2 * sqrt( E * b^3*h/12 * E/2/(1+v) * b^3*h/3 );
        b1 = 2.78088772399498*2/Le^2;

        theta2= asin(1/gamma);
        vB = L1*gamma*(deg2rad(theta)+asin(1/gamma))  %a way to increase A1

        syms A1 z
        phi = sqrt(Le-z).*A1.*besselj(1/4,b1/2.*(Le-z)^2);
        U = Pcr2^2/2/E/ (b^3*h/12) *int( sin(phi).^2.*(Le-z).^2 ,z,0,Le) + E/2/(1+v)*b^3*h/3 /2*int( diff(phi,z,1).^2 ,z,0,Le);%sin(phi) to increase psi?

        T = Pcr2*vB;

        eq1 = U == T;
        A1 = double(vpasolve(eq1, A1, 0.100332869113339))  % print A1 value
%         du_dz0(i,2)=A1;
        %4
        Pcr2 = 5.5618/L^2 * sqrt( E * b^3*h/12 * E/2/(1+v) * b^3*h/3 );%

        phi =@(z) sqrt(Le-z).*A1.*besselj(1/4,b1/2.*(Le-z)^2);%

        % f = fplot(phi,[0 Le]);
        % x1 = f.XData;
        % y1 = f.YData;
        %8
        syms z
        phi1 = sqrt(Le-z).*A1.*besselj(1/4,b1/2.*(Le-z)^2);
        du_dz = @(a) Pcr2/E/(b^3*h/12) *int(phi1*(Le-z),z, 0, a);% * cos(subs(phi,z,a)); what is the point of adding cos(phi)

        du_dz0(i,j) = double(du_dz(Le));
    end
end
for i = 1: length(theta0)
    for j = 1:length(gamma0)
        theta=theta0(i);
        gamma = gamma0(j);

        Ebarr0(i,j)=shellBarrEner(theta,gamma);
    end
end
%% 21. Sep 6 2022 Calculate the biped robot psi_l and U_bar, borrowing from 20
clear

l0 = [129.1];%l0 in mm
D0 = [16];

v = 0.34;%0.43; 

h = 15;%mm
E = 1738;%1738; % mm, N, MPa, s
b = 0.381; %mm Pcr2 calculated from above-above section 0.1227


du_dz0 = zeros(length(l0),length(D0));
Ebarr0 = zeros(length(l0),length(D0));

for i = 1: length(l0)
    for j = 1:length(D0)
        
        l=l0(i);
        D=D0(j);

        L = l;
        Le = l;
        Pcr2 = 2.78088772399498*2/L^2 * sqrt( E * b^3*h/12 * E/2/(1+v) * b^3*h/3 );
        b1 = 2.78088772399498*2/L^2;

        vB = D/2;  %a way to increase A1,

        syms A1 z
        phi = sqrt(Le-z).*A1.*besselj(1/4,b1/2.*(Le-z)^2);
        U = Pcr2^2/2/E/ (b^3*h/12) *int( sin(phi).^2.*(Le-z).^2 ,z,0,Le) + E/2/(1+v)*b^3*h/3 /2*int( diff(phi,z,1).^2 ,z,0,Le);%sin(phi) to increase psi?

        T = Pcr2*vB;

        eq1 = U == T;
        A1 = double(vpasolve(eq1, A1, 0.100332869113339))  % print A1 value
%         du_dz0(i,2)=A1;
        %4
        Pcr2 = 5.5618/L^2 * sqrt( E * b^3*h/12 * E/2/(1+v) * b^3*h/3 );%

        phi =@(z) sqrt(Le-z).*A1.*besselj(1/4,b1/2.*(Le-z)^2);%

        % f = fplot(phi,[0 Le]);
        % x1 = f.XData;
        % y1 = f.YData;
        %8
        syms z
        phi1 = sqrt(Le-z).*A1.*besselj(1/4,b1/2.*(Le-z)^2);
        du_dz = @(a) Pcr2/E/(b^3*h/12) *int(phi1*(Le-z),z, 0, a);% * cos(subs(phi,z,a)); what is the point of adding cos(phi)

        du_dz0(i,j) = double(du_dz(Le));
    end
end
for i = 1: length(l0)
    for j = 1:length(D0)
        l=l0(i);
        D = D0(j);

        Ebarr0(i,j)=shellBarrEner_LD(l,D);

    end
end
%% 22. Calculating lateral displacement or deflection from du_dz, borrowing from 20
clear

l0 = [187.49/2];%l0 in mm, half ribbon length
D0 = [20];

%gamma0 = [4];%(2 5)=53.34, (2 6)=57.49, (0 6)=57.49, (0 4 20)=48.84, (0 4 20)=48.84, (0 4 20 0.43)=48.09,
             %(0 4 15)=42.1  %(0 4 0)=8e-30

v = 0.43;%0.43 decrease the value a bit 

h = 15;%mm
E = 1738;%1738; % mm, N, MPa, s
b = 0.381; %mm Pcr2 calculated from above-above section 0.1227. 0.381 mm = 0.015 inch


du_dz0 = zeros(length(l0),length(D0));
Ebarr0 = zeros(length(l0),length(D0));

u0 = zeros(length(l0),length(D0));

for i = 1: length(l0)
    for j = 1:length(D0)

        l=l0(i);
        D=D0(j);

        L = l;
        Le = l;
        Pcr2 = 2.78088772399498*2/L^2 * sqrt(E * b^3*h/12 * E/2/(1+v) * b^3*h/3 );
        b1 = 2.78088772399498*2/Le^2;

        %vB = L1*gamma*(deg2rad(theta)+asin(1/gamma));  %a way to increase A1, this is the displacement, a curved path
        vB = D/2; %remember it is a half of the D

        syms A1 z
        phi = sqrt(Le-z).*A1.*besselj(1/4,b1/2.*(Le-z)^2);
        U = Pcr2^2/2/E/ (b^3*h/12) *int( sin(phi).^2.*(Le-z).^2 ,z,0,Le) + E/2/(1+v)*b^3*h/3 /2*int( diff(phi,z,1).^2 ,z,0,Le);%sin(phi) to increase psi?

        T = Pcr2*vB;

        eq1 = U == T;
        A1 = double(vpasolve(eq1, A1, 0.100332869113339))  % print A1 value
%         du_dz0(i,2)=A1;
        %4
        Pcr2 = 5.5618/L^2 * sqrt( E * b^3*h/12 * E/2/(1+v) * b^3*h/3 );%

        phi =@(z) sqrt(Le-z).*A1.*besselj(1/4,b1/2.*(Le-z)^2);%

        % f = fplot(phi,[0 Le]);
        % x1 = f.XData;
        % y1 = f.YData;
        %8
        syms z
        phi1 = sqrt(Le-z).*A1.*besselj(1/4,b1/2.*(Le-z)^2);
        du_dz = @(a) Pcr2/E/(b^3*h/12) *int(phi1*(Le-z),z, 0, a);% * cos(subs(phi,z,a)); what is the point of adding cos(phi)
        
        syms aa
        u =@(bb) int( du_dz(aa), aa, 0, bb);

        du_dz0(i,j) = double(du_dz(Le));
        u0(i,j) = double(u(Le));
    end
end

%% 23. Mar 7 2023 calcalating ornithopter clap-fling angle, and required forces, borrow from 22
clear

l0 = 87.5%[150];%l0 in mm, l0 is the "half" ribbon length
D0 = 17.145%[25];
L1 = 20%[44]; % 30 or 12.5+7.5 

v = 0.38;%0.26;%0.43 decrease the value a bit, 0.26 is for CFRP

h = 15;%mm  angle is irrelevant to h,E,b
E = 1738;%1738; % mm, N, MPa, s
b = 0.381; %mm Pcr2 calculated from above-above section 0.1227. 0.381 mm = 0.015 inch


du_dz0 = zeros(length(l0),length(D0));
du_dz_L1 = zeros(length(l0),length(D0));
Ebarr0 = zeros(length(l0),length(D0));

u0 = zeros(length(l0),length(D0));
u_L1 = zeros(length(l0),length(D0));

for i = 1: length(l0)
    for j = 1:length(D0)

        l=l0(i);
        D=D0(j);

        L = l;
        Le = l;
        Pcr2 = 2.78088772399498*2/L^2 * sqrt(E * b^3*h/12 * E/2/(1+v) * b^3*h/3 );
        b1 = 2.78088772399498*2/Le^2;

        %vB = L1*gamma*(deg2rad(theta)+asin(1/gamma));  %a way to increase A1, this is the displacement, a curved path
        vB = D/2; %remember it is a half of the D

        syms A1 z
        phi = sqrt(Le-z).*A1.*besselj(1/4,b1/2.*(Le-z)^2);
        U = Pcr2^2/2/E/ (b^3*h/12) *int( sin(phi).^2.*(Le-z).^2 ,z,0,Le) + E/2/(1+v)*b^3*h/3 /2*int( diff(phi,z,1).^2 ,z,0,Le);%sin(phi) to increase psi?

        T = Pcr2*vB;

        eq1 = U == T;
        A1 = double(vpasolve(eq1, A1, 0.100332869113339))  % print A1 value
%         du_dz0(i,2)=A1;
        %4
        Pcr2 = 5.5618/L^2 * sqrt( E * b^3*h/12 * E/2/(1+v) * b^3*h/3 );%

        phi =@(z) sqrt(Le-z).*A1.*besselj(1/4,b1/2.*(Le-z)^2);%
        %8  % to calculate overall displacement
        syms z
        phi1 = sqrt(Le-z).*A1.*besselj(1/4,b1/2.*(Le-z)^2);
        du_dz = @(a) Pcr2/E/(b^3*h/12) *int(phi1*(Le-z),z, 0, a);% * cos(subs(phi,z,a)); what is the point of adding cos(phi)
        
        syms aa
        u =@(bb) int( du_dz(aa), aa, 0, bb);

        du_dz0(i,j) = double(du_dz(Le));
        u0(i,j) = double(u(Le));

        syms z   %to calculate L1 point displacement
        phi1 = sqrt(Le-z).*A1.*besselj(1/4,b1/2.*(Le-z)^2);
        du_dz = @(a) Pcr2/E/(b^3*h/12) *int(phi1*(Le-z),z, 0, a);% * cos(subs(phi,z,a)); what is the point of adding cos(phi)
        
        syms aa
        u =@(bb) int( du_dz(aa), aa, 0, bb);

        du_dz_L1(i,j) = double(du_dz(L1));
        u_L1(i,j) = double(u(L1));

    end
end

for i = 1: length(l0)
    for j = 1:length(D0)
        l=l0(i);
        D = D0(j);

        Ebarr0(i,j)=shellBarrEner_LD(l,D);

    end
end

%% 24. Apri 8 2023 calcalating ornithopter <Fish 0.79CFRP + B24CLM> angle, and required forces, borrow from 23
clear

l0 = [137];%l0 in mm, l0 is the "half" ribbon length
D0 = [10];
L1 = [44]; % 

v = 0.26;%0.26;%0.43 decrease the value a bit, 0.26 is for CFRP

h = 10;%mm  angle is irrelevant to h,E,b
E = 64000;%1738; % mm, N, MPa, s
b = 0.5 ;%0.79; %mm Pcr2 calculated from above-above section 0.1227. 0.381 mm = 0.015 inch


du_dz0 = zeros(length(l0),length(D0));
du_dz_L1 = zeros(length(l0),length(D0));
Ebarr0 = zeros(length(l0),length(D0));

u0 = zeros(length(l0),length(D0));
u_L1 = zeros(length(l0),length(D0));

for i = 1: length(l0)
    for j = 1:length(D0)

        l=l0(i);
        D=D0(j);

        L = l;
        Le = l;
        Pcr2 = 2.78088772399498*2/L^2 * sqrt(E * b^3*h/12 * E/2/(1+v) * b^3*h/3 );
        b1 = 2.78088772399498*2/Le^2;

        %vB = L1*gamma*(deg2rad(theta)+asin(1/gamma));  %a way to increase A1, this is the displacement, a curved path
        vB = D/2; %remember it is a half of the D

        syms A1 z
        phi = sqrt(Le-z).*A1.*besselj(1/4,b1/2.*(Le-z)^2);
        U = Pcr2^2/2/E/ (b^3*h/12) *int( sin(phi).^2.*(Le-z).^2 ,z,0,Le) + E/2/(1+v)*b^3*h/3 /2*int( diff(phi,z,1).^2 ,z,0,Le);%sin(phi) to increase psi?

        T = Pcr2*vB;

        eq1 = U == T;
        A1 = double(vpasolve(eq1, A1, 0.100332869113339))  % print A1 value
%         du_dz0(i,2)=A1;
        %4
        Pcr2 = 5.5618/L^2 * sqrt( E * b^3*h/12 * E/2/(1+v) * b^3*h/3 );%

        phi =@(z) sqrt(Le-z).*A1.*besselj(1/4,b1/2.*(Le-z)^2);%
        %8  % to calculate overall displacement
        syms z
        phi1 = sqrt(Le-z).*A1.*besselj(1/4,b1/2.*(Le-z)^2);
        du_dz = @(a) Pcr2/E/(b^3*h/12) *int(phi1*(Le-z),z, 0, a);% * cos(subs(phi,z,a)); what is the point of adding cos(phi)
        
        syms aa
        u =@(bb) int( du_dz(aa), aa, 0, bb);

        du_dz0(i,j) = double(du_dz(Le));
        u0(i,j) = double(u(Le));

        syms z   %to calculate L1 point displacement
        phi1 = sqrt(Le-z).*A1.*besselj(1/4,b1/2.*(Le-z)^2);
        du_dz = @(a) Pcr2/E/(b^3*h/12) *int(phi1*(Le-z),z, 0, a);% * cos(subs(phi,z,a)); what is the point of adding cos(phi)
        
        syms aa
        u =@(bb) int( du_dz(aa), aa, 0, bb);

        du_dz_L1(i,j) = double(du_dz(L1));
        u_L1(i,j) = double(u(L1));

    end
end

for i = 1: length(l0)
    for j = 1:length(D0)
        l=l0(i);
        D = D0(j);

        Ebarr0(i,j)=shellBarrEner_LD(l,D);

    end
end
%% #25 plot besselJ Aug 23
x = 0:0.01:10;
y = besselj(1/4,x);
plot(x,y)
%% #26 Aug 23, Limit Cycles, nonlinear Control Systems
clear all
%program to analyse Van der Pol oscillator
mu = 1;
[x1, x2] = meshgrid(-4:1/2/2:4, -4:1/2/2:4);
dx1 =x2;
dx2 = -mu*(x1.^2 -1) .*x2 -x1; % this is van der pol
% dx2 = -mu*x2 -x1; % this is not
r = sqrt(dx1.^2 + dx2.^2);
quiver(x1, x2, dx1./r, dx2./r, 1/2, 'LineWidth', 1);

hold on;
axis equal;
set(gca, 'Fontsize', 20);
axis( [-4 4 -4 4] );

while true
    x0 = ginput(1);
    tspan = [0 60];
    [t,x] = ode45(@(t,x) odefcn(x, mu), tspan, x0);
    plot(x(:, 1), x(:,2),'LineWidth',3)
end




%% #27 Aug 23 Solving improved CPG model from a paper of ？
% fplot(t alpha, [0 0.78], 'b')

% ka = 10;
% kb = 7;
% km = 7;
% 
% tau = 50;
% 
% % sym t
% % u = A .* sin(omega .* t)
% 
% [t,y] = ode45(@vdp1,[0 20],[2; 0]);
% 
% 
% plot(t,y(:,1),'-o',t,y(:,2),'-o')
% title('Solution of van der Pol Equation (\mu = 1) with ODE45');
% xlabel('Time t');
% ylabel('Solution y');
% legend('y_1','y_2')
%% #28 Aug 23 plotting smoothed curve for ICRA 2024 Paper
clear all

f = 1/0.760; %Hz , unit is s?!!
omega = f *2 *pi();

R = 0.5;
B = 2;
M= 57.09; % bending angle in deg

t = linspace(0.325,3.375,367);% 0.325 + 4T = 3.365

% alpha = @(t) M *  sin( omega * t);

alpha = @(t) M * tanh(B* sin( omega * t)) / tanh(B);

v1 = alpha(t);

alpha2 = @(x) 56.87683*sin(pi*(x-0.01243483)/0.380);

v2 = alpha2(t);
plot(t,v2,'o')
%% #29 Sept 23 Fitting HCM undulation to VDP oscillation
open_system('vdp');
set_param('vdp/Mu','Gain','1')
sim('vdp');
open_system('vdp/Scope');
%% #30 Oct 23 Calculating for the RoboSoft 2024
% calculate F_max <- Energy , f_m, HCM <- directly from the equ of t*
clear

l0 = 137;%[67  87  107  127  147 167 ];%l0 in mm, l0 is the "half" ribbon length
D0 = 10%[7.5 10 12.5 15 17.5 20 22.5]; %
% calculation of u(L1) is not needed, just assume the rotation angle = 1 rad

v0 = [0.26 0.26 0.3 0.301  0.301 0.3];%0.26;%0.43 decrease the value a bit, 0.26 is for CFRP
h0 = [10 15 10 15  15 10];%mm  angle is irrelevant to h,E,b
E0 = [64e3 64e3 200e3 1738 1738 200e3];%1738; % mm, N, MPa, s 64000 for CFRP
b0 = [0.5 0.5 0.15 0.762  0.381 0.5];%79; %mm Pcr2 calculated from above-above section 0.1227. 0.381 mm = 0.015 inch
rho0 = [1.6e-9 1.6e-9 7.8e-9 1.25e-9 1.25e-9 7.8e-9]; % unit in t/mm3
horn_L0 = [25 25 25 10 10 25]; %mm
L10 = [44 44 44 30 12.5 44];

Ebarr0 = zeros(length(v0),length(l0),length(D0)); %unit in mJ
Tor_m_HCM0 = zeros(length(v0),length(l0),length(D0));

TimeStar0 = zeros(length(v0),length(l0), 1);
f_m_HCM0 = zeros(length(v0),length(l0), 1);

du_dz0 = zeros(length(v0),length(l0),length(D0));
du_dz_L1 = zeros(length(v0),length(l0),length(D0));

u0 = zeros(length(v0),length(l0),length(D0));
u_L1 = zeros(length(v0),length(l0),length(D0));

for m = 6:6% length(v0)
    
    v = v0(m);
    h = h0(m);
    E = E0(m);
    b = b0(m);
    rho = rho0(m);
    horn_L = horn_L0(m);
    L1 = L10(m);
    
    for i = 1: length(l0)
        for j = 1:length(D0)            
            
            l=l0(i);
            D=D0(j);

            L = l;
            Le = l;
            Pcr2 = 2.78088772399498*2/L^2 * sqrt(E * b^3*h/12 * E/2/(1+v) * b^3*h/3 );
            b1 = 2.78088772399498*2/Le^2;

            %vB = L1*gamma*(deg2rad(theta)+asin(1/gamma));  %a way to increase A1, this is the displacement, a curved path
            vB = D/2; %remember it is a half of the D

            syms A1 z
            phi = sqrt(Le-z).*A1.*besselj(1/4,b1/2.*(Le-z)^2);
            U = Pcr2^2/2/E/ (b^3*h/12) *int( sin(phi).^2.*(Le-z).^2 ,z,0,Le) + E/2/(1+v)*b^3*h/3 /2*int( diff(phi,z,1).^2 ,z,0,Le);%sin(phi) to increase psi?

            T = Pcr2*vB;

            eq1 = U == T;
            A1 = double(vpasolve(eq1, A1, 0.100332869113339))  % print A1 value
    %         du_dz0(i,2)=A1;

            Pcr2 = 5.5618/L^2 * sqrt( E * b^3*h/12 * E/2/(1+v) * b^3*h/3 );%

            phi =@(z) sqrt(Le-z).*A1.*besselj(1/4,b1/2.*(Le-z)^2);%

            %8  % to calculate overall displacement
            syms z
            phi1 = sqrt(Le-z).*A1.*besselj(1/4,b1/2.*(Le-z)^2);
            du_dz = @(a) Pcr2/E/(b^3*h/12) *int(phi1*(Le-z),z, 0, a);% * cos(subs(phi,z,a)); what is the point of adding cos(phi)

            syms aa
            u =@(bb) int( du_dz(aa), aa, 0, bb);

            du_dz0(m,i,j) = double(du_dz(Le));
            u0(m,i,j) = double(u(Le));

            syms z   %to calculate L1 point displacement
            phi1 = sqrt(Le-z).*A1.*besselj(1/4,b1/2.*(Le-z)^2);
            du_dz = @(a) Pcr2/E/(b^3*h/12) *int(phi1*(Le-z),z, 0, a);% * cos(subs(phi,z,a)); what is the point of adding cos(phi)

            syms aa
            u =@(bb) int( du_dz(aa), aa, 0, bb);

            du_dz_L1(m,i,j) = double(du_dz(L1));
            u_L1(m,i,j) = double(u(L1));

        end

        TimeStar0(m,i) = (2*l)^2/b/sqrt(E/rho); %unit in s, calculating the timescale of each half length
        f_m_HCM0(m,i) = 1 ./ (TimeStar0(i) .* 2); %unit in s 

    end

    for i = 1: length(l0)
        for j = 1:length(D0)
            l=l0(i);
            D = D0(j);

            Ebarr0(m,i,j)=shellBarrEner_LD(l,D);

            Tor_m_HCM0(m,i,j) = Ebarr0(m,i,j) * 2 / u_L1(m,i,j) * horn_L ; % unit in N*mm, not accurate??? like

        end
    end
end
%3300 to 150, 20 times difference between needed and offered by the servo

%% #30 similar to #29, calculating the frequency for l and h, RoboSoft2024
clear

l0 = [67 77 87 97 107 117 127 137 147 157 167];%l0 in mm, l0 is the "half" ribbon length
% calculation of u(L1) is not needed, just assume the rotation angle = 1 rad

v = 0.3;%0.26;%0.43 decrease the value a bit, 0.26 is for CFRP
E = 1738;%1738; % mm, N, MPa, s
b = 0.762;%79; %mm Pcr2 calculated from above-above section 0.1227. 0.381 mm = 0.015 inch
rho = 1.25e-9; % unit in t/mm3

TimeStar0 = zeros(length(l0), 1);
f_m_HCM0 = zeros(length(l0), 1);

for i = 1: length(l0)

    l=l0(i);
    TimeStar0(i) = (2*l)^2/b/sqrt(E/rho); %unit in s, calculating the timescale of each half length
    f_m_HCM0(i) = 1 ./ (TimeStar0(i) .* 2); %unit in s 
    
end                        %%%%%% draw the surface so that the points are on them,

%% Function library
function dydt = vdp1(t,y)  %function only damp y1 y2??? yes 
kb = 1;
B = 0;
%VDP1  Evaluate the van der Pol ODEs for mu = 1
%
%   See also ODE113, ODE23, ODE45.

%   Jacek Kierzenka and Lawrence F. Shampine
%   Copyright 1984-2014 The MathWorks, Inc.

dydt = [y(2); kb * (0.25*kb*(B-y(1))-y(2)) ];
end 

function dxdt = odefcn(x,mu)
dxdt = zeros(2,1);
dxdt(1) = x(2);
dxdt(2) = -mu*(x(1)^2 - 1)*x(2) - x(1);
end

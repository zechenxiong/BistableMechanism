function [Ebarr,Ebarr_, Emax_,Emax]= shellBarrEner(theta, gamma, v, L1, b, h, E ) %theta in deg, all else in mm and MPa
%function [U_s_,D_U_,D_U_2,U_m_,D_U_3]= shellBarrEner(theta, g_s ), h =15
%b=0.381
%==
% g_s = L2/L1  shape factor, or use gamma
% D_U_ = D_U*L1/E/I_eta, U bar, normalized energy barrier
% in the FEM, L1 = 12.5, E=2000, Ieta=15*0.381^3/12,
% DU_=DU*12.5/2000/15/0.381^3*12=DU*0.0904055468

%L1 = 12.5*1;
L = L1*(gamma+1);
%b = 0.381;
%h = 15;
% v = 0.34; %used to be 0.34, now 0.4

% E = 1730; %MPa PETG modulus, Steel modulus 200GPa
G = E/2/(1+v);
J = b^3*h/3; 
Ieta = b^3*h/12;
C = G*J; % torsional rigidity

Pcr = 2*2.78088772399498/L^2 * sqrt( E * Ieta * G * J );
vB = L1*gamma*(deg2rad(theta)+asin(1/gamma));

Es = 2*Pcr*vB;% elastic energy stored after prestressing of a HCM
Es_= Es*L1/E/h/b^3*12;

Emax = Es*4;%how to prove??? then E_barr is 7 times the Es
Emax_ = Emax*L1/E/h/b^3*12;

Ebarr = Emax/4*3;  %unit is mJ
Ebarr_ = Emax_/4*3;


%% prove, according to the paper "configuration evolution...", the symmetric snapping can have 
% y|x=l = 0 as the unstable bifurcation point, which correspond to the E_max state, when this happen,
% the amplitude of the second-order buckling profile is half of the
% amplitude of the first order buckling profile, ie, the initial profile.
% thus the elastic energy
% prestressed profile: y= A cos(2·pi/l·x)
% turns out it should be 4 times and 3 times????
% how to prove the new amplitude is A/sqrt(2)???



% D_U_ = 2*5.5618... %normalized snapping energy barrier
%        * ( (1+g_s-h/4/L1)^2/(g_s-3*h/4/L1)^2 -1)...
%        * sqrt( 2/(1+v) )...
%        * ( 1/cos(theta)*1/(1+g_s-h/4/L1) + g_s/(1+g_s-h/4/L1)*tan(theta) )...
%        / (1+g_s-h/4/L1);
%    
% U_s_ = 2*5.5618...#normalized stored elastic energy
%        * ( 1/(1+g_s-h/4/L1) )^2 ...
%        * sqrt( 2/(1+v) )...
%        * (1+g_s*sin(theta))/cos(theta);
%   
% D_U_2 = 2*5.5618...
%         * ( 1/(g_s-3*h/4/L1)^2 - 1/(1+g_s-h/4/L1)^2 )...
%         * sqrt(2/(1+v))...
%         * (1+g_s*sin(theta))/cos(theta);
%     
% U_m_ = 2*5.5618... 
%        * 1/(g_s-3*h/4/L1)^2 ...
%        * sqrt( 2/(1+v) )...
%        * (1+g_s*sin(theta))/cos(theta);
% 
% D_U_3 = U_m_-U_s_;%another snapping energy barrier?

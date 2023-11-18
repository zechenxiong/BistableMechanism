function [U_barr]= shellBarrEner_LD(l, D ) %theta in deg
%function [U_s_,D_U_,D_U_2,U_m_,D_U_3]= shellBarrEner(theta, g_s )
% g_s = L2/L1  shape factor, or use gamma
% D_U_ = D_U*L1/E/I_eta, U bar, normalized energy barrier
% in the FEM, L1 = 12.5, E=2000, Ieta=15*0.381^3/12,
% DU_=DU*12.5/2000/15/0.381^3*12=DU*0.0904055468
% what is the fucking Emax, why is it 8 times the Es? changed to 4
% times.---no 8 times is correct, because the HCM has two props, so 4*2
% what is the unit of the U_barrr?

v = 0.30; %0.26 is CFRP
h = 15;
E = 200e3;%64000;%1738; % in MPa
b = 0.15;%0.381;


L = l; % in mm 
G = E/2/(1+v);
J = b^3*h/3; 
Ieta = b^3*h/12; % in mm^4
C = G*J; % torsional rigidity

Pcr = 2*2.78088772399498/L^2 * sqrt( E * Ieta * G * J );
vB = D;

Es = Pcr*vB;% elastic energy stored after prestressing, unit in N*mm

U_barr = 3* Es; % in mN, mJ? This is for the whole
% Emax = Es*4;
% Emax_ = Emax*L1/E/h/b^3*12;



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

clear all

pwr=[300 500]; %[pwr]=mega-volt-amperes
weight=[225 325]*0.9071847; %[weight]=mT
p=polyfit(pwr,weight,1);
weight=p(1)*(2.4420e+21)/(1e6)+p(2)

% 20kV capacitors
% weight=p(1)*(4.9038e+20)/(1e6)+p(2)
% tot_cap=110e-3; %[tot_cap]=F
% cap=6.8e-9; %[cap]=F
% num_caps=tot_cap/cap
% capD=1.4; %[capD]=in
% capT=0.65; %[capT]=in
% vol_cap=(pi/4)*1.4^2*0.65;
% rho=72; %[rho]=g/in^3
% cap_mass=rho*vol_cap*num_caps
% % mass in g
% cap_mass*1e-6
% % mass in mT

% 10kV capacitors
tot_cap=450e-3; %[tot_cap]=F
cap=562e-12; %[cap]=F
num_caps=tot_cap/cap
capL=1.3*0.03937008; %[capD]=in
capW=0.6*0.03937008; %[capD]=in
capT=0.222*0.03937008; %[capT]=in
vol_cap=capL*capW*capT;
rho=67; %[rho]=g/in^3
cap_mass=rho*vol_cap*num_caps
% mass in g
cap_mass*1e-6
% mass in mT
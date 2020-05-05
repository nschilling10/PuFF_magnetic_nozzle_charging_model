%% Nathan Schilling
% Transformer circuit model wrapper script
% Simply looks at specifed values of L1 & L_ratio. No trade study
% Seed current is 2.2e5 A to get correct circuit energy & capacitor energy
% 02/19/19
clear all
close all

test.L1=1e-6;
test.Lratio=1;
test.graphDisplay=true;
[E_gain,E_circ] = circuitModelFunv2(test)

%% Impacts on thrust & Isp
Plasma_thermal_energy_calc
E_avail=E_therm-E_circ
eta=0.2;
g_0=9.81;
f=10;
V_e=sqrt(2*eta*E_avail/(m_T*1e-3));
Isp=V_e/g_0
T=f*V_e*m_T

E_circ=0;
E_avail=E_therm-E_circ;
V_e=sqrt(2*eta*E_avail/(m_T*1e-3));
Isp=V_e/g_0
T=f*V_e*m_T
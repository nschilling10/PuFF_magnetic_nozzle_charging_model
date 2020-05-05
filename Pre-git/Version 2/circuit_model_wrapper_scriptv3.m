%% Nathan Schilling
% Takes optimized results from trade study and generates output file
% Also determines results on specific impulse and thrust
% 02/19/19
clear all
close all
format long

test.L1=20e-6;
test.L2=1e-6;
test.L0=130e-6;
test.I0=4e5;

test.R1=10;
test.R2=10;

test.graphDisplay=true;
[E_gain,E_circ] = circuitModelFunv4(test)
Gain=E_gain/E_circ
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
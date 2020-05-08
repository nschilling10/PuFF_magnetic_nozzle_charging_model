%% Model Wrapper script
% Does a run of the circuit model for a single set of values
% Not the parametric study code
% Once you do a parametric study, run the code again here to see the
% results
% Nathan Schilling
% 05/05/2020

clear all
close all

lectureNotes=false; 
PUFF=false;
reactorGradeCaps=false;
%% Physical Constants
g0=9.80665; %[g0]=m/s^2
Ru=8.31446261815324e3; %[Ru]=J/k-mol-K
k_b=8.617333262145e-5; %[k_b]=eV/K
e=1.602176634e-19; %[e]=J/eV
eps0=8.8541878128e-12; %[eps0]=SI units
mu_0=4*pi*1e-7; %[mu_0]=SI units
h=6.62607015e-34; %[h]=J*s
m_a=1.6605e-27; %[m_a]=kg/u
m_e=9.10938356e-31; %[m_e]=kg
%% Unit conversion
eV = 1/k_b;
keV = 1000*eV;
%% inputs
plasmaInps.m_propellant = .01;
if lectureNotes %case used in lecture notes, 40 kV charged to caps
    primary circuit (1) inputs
    circInps.L1 = 100e-9;
    circInps.R1 = 5e-3;
    
    %secondary circuit (2) inputs including reactor capacitance C
    circInps.L2 = 100e-9;
    circInps.R2 = 5e-3;
    circInps.C = 100e-6;
    
    %plasma inputs
    plasmaInps.MW = 4;
    Rgas = Ru/MW;
    plasmaInps.T0 = .5; %[T0]=keV
    plasmaInps.g = 1.3;
    
    %transformer inputs
    circInps.mu_r = 1;
    circInps.k=0.9;
    
    circInps.rT = .1;
    circInps.NT1 = 10;
    circInps.NT2 = NT1/10;
    circInps.lT1 = rT*1;
    circInps.lT2 = rT*3;
    
    %flux compression generator inputs
    circInps.Nfcc = 10;
    circInps.Rfcc = .95; %1.1;  %coil radius
    
    %-----------initial conditions
    circInps.I10 = .51e6; %Amps, initial current in primary
    circInps.I20 = 0;
%     circInps.V10 = 0;
    circInps.V20 = 0;  %voltage on reactor caps
elseif PUFF  %PUFF
    plasmaInps.m_propellant = .125;
    %primary circuit (1) inputs
    circInps.L1 = 100e-9;
    circInps.R1 = 0*1e-5 ;
    
    %secondary circuit (2) inputs including reactor capacitance C
    circInps.L2 = 10e-9;
    circInps.R2 = 5e-3 * 0;
    circInps.C = .01;  %because we need 50 MJ or more
    
    %plasma inputs
    plasmaInps.MW = 4;
    Rgas = Ru/plasmaInps.MW;
    plasmaInps.T0 = 1*.98; %[T0]=keV
    plasmaInps.g = 1.3;
    
    %transformer inputs
    circInps.mu_r = 1;
    circInps.k=0.9;
    
    circInps.rT = .1;
    circInps.NT1 = 30;
    circInps.NT2 = 1.7; %NT1/10;
    circInps.lT1 = circInps.rT*1;
    circInps.lT2 = circInps.rT*4;
    
    %flux compression generator inputs
    circInps.Nfcc = 9;
    circInps.Rfcc = 18;
    
    %-----------initial conditions
    circInps.I10 = .75e6; %Amps, initial current in primary
    circInps.I20 = 0;
    circInps.V10 = 0;
    circInps.V20 = 0;  %voltage on reactor caps
elseif reactorGradeCaps  %reactor grade caps
    
    plasmaInps.m_propellant = .01;
    %primary circuit (1) inputs
    circInps.L1 = 100e-9;
    circInps.R1 = 5e-3 ;
    
    %secondary circuit (2) inputs including reactor capacitance C
    circInps.L2 = 100e-9;
    circInps.R2 = 5e-3 * 1;
    circInps.C = 400e-6;
    
    %plasma inputs
    plasmaInps.MW = 4;
    plasmaInps.T0 = 1; %[T0]=keV
    plasmaInps.g = 1.3;
    
    %transformer inputs
    circInps.mu_r = 1;
    
    circInps.rT = .1;
    AT = rT^2 * pi;
    circInps.NT1 = 30;
    circInps.NT2 = 1.7; %NT1/10;
    circInps.lT1 = circInps.rT*1;
    circInps.lT2 = circInps.rT*4;
    
    %flux compression generator inputs
    circInps.Nfcc = 10;
    circInps.Rfcc = 3;
    
    %-----------initial conditions
    circInps.I10 = .77e6; %Amps, initial current in primary
    circInps.I20 = 0;
%     circInps.V10 = 0;
    circInps.V20 = 0;  %voltage on reactor caps
else % your stuff goes here
    
    plasmaInps.m_propellant = .150;
    %primary circuit (1) inputs
    circInps.L1 = 100e-9;
    circInps.R1 = 1e-5 * 0;
    
    %secondary circuit (2) inputs including reactor capacitance C
    circInps.L2 = 10e-9;
    circInps.R2 = 5e-3 * 0;
    circInps.C = .01;  %because we need 50 MJ or more
    
    %plasma inputs
    plasmaInps.MW = 235;
    Rgas = Ru/plasmaInps.MW;
    plasmaInps.T0 = 50;
    plasmaInps.g = 1.3;
    
    %transformer inputs
    circInps.mu_r = 1;
    circInps.k=0.9;
    
    circInps.rT = .1;
    circInps.NT1 = 300;
    circInps.NT2 = 1; %NT1/10;
    circInps.lT1 = circInps.rT*1;
    circInps.lT2 = circInps.rT*1;
    
    %flux compression generator inputs
    circInps.Nfcc = 20;
    circInps.Rfcc = 17;
    
    %-----------initial conditions
    circInps.I10 = 0.77e6; %Amps, initial current in primary
    circInps.I20 = 0;
%     circInps.V10 = 0;
    circInps.V20 = 0;  %voltage on reactor caps
end
%%
graphDisplay=true;
chargingModelOutputs=charging_nozzle_model(circInps,plasmaInps,graphDisplay);
%% Nozzle parameters
Vf=-chargingModelOutputs.V_plasma(end);
Isp=Vf/g0
Impulse_bit=plasmaInps.m_propellant*Vf
%% Radiation loads
Rgas = Ru/plasmaInps.MW;
Yield=3*plasmaInps.m_propellant*Rgas*plasmaInps.T0*keV %[Yield]=J
E_fiss=207*1e6*e; %[E_fiss]=J/fiss
num_fiss=Yield/E_fiss
neut_per_fiss=2.43;
tot_neut=neut_per_fiss*num_fiss;
neturon_fluence_per_steradian=tot_neut/(4*pi)

neut_Energy_per_fiss=5*1e6*e;
tot_neut_Energy=neut_Energy_per_fiss*num_fiss;
neturon_Energy_fluence_per_steradian=tot_neut_Energy/(4*pi)

photonEnergy_per_fiss=14*1e6*e;
tot_photon_energy=photonEnergy_per_fiss*num_fiss;
tot_photon_energy_per_steradian=tot_photon_energy/(4*pi)
function [eta_KE,eta_circ,eta_sys]=efficiencyCalcFun(energies,plasmaInps)
%% Efficiency Calculation function
% Nathan Schilling
% 05/05/2020
% Calculates kinetic energy, circuit, and system efficiencies for a given
% nozzle/circuit design. 

eV=11605;
Ru=8314.5;
Rgas=Ru/plasmaInps.MW; %this requires temeprature in KELVIN

E_plasma=energies.E_plasma;
E_therm0=3*plasmaInps.m_propellant*Rgas*plasmaInps.T0*1e3*eV*1e-6;
E_cap=energies.E_cap;
E_ohmic_losses=energies.E_ohmic_losses;

% ---------- Efficiency Calculation
eta_KE=E_plasma(end)/(E_therm0-E_cap(end));

eta_circ=1-E_ohmic_losses(end)/(E_ohmic_losses(end)+E_cap(end));

eta_sys=(E_plasma(end)+E_cap(end))/E_therm0;

end
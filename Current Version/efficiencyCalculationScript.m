function [eta_KE,eta_circ,eta_sys]=efficencyCalcFunction(energies)

%% Efficiency Calculation function
% Nathan Schilling
% 05/05/2020
% Calculates system, circuit, and kinetic energy efficiencies for a given
% nozzle/circuit design. Requires energy calculation script to be in the
% same working directory.

E_plasma=energies.E_plasma;
E_therm0=energies.E_therm(1);
E_cap=energies.E_cap;
E_ohmic_losses=energies.E_ohmic_losses;

% ---------- Efficiency Calculation

eta_KE=E_plasma(end)/(E_therm0-E_cap(end));

eta_circ=1-E_ohmic_losses(end)/(E_ohmic_losses(end)+E_cap(end));

eta_sys=(E_plasma(end)+E_cap(end))/E_therm0;

end
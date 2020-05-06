function energies=calcCompEnergiesFun(modelOutputs, circInps, plasmaInps)
%% Energy Calculation function
% Nathan Schilling
% 05/05/2020
% Holds all the necessary ancilary calcualtions needed to calculate various
% energies in the circuit (Fcc energy, plasma KE, load capacitor PE,
% resisitve losses, thermal energy using isentropic flow assumption)
% INPUTS
% modelOutputs is assumed to come from the charging circuit model - it
% should have t, I1, I2, Vfcc, Vcap, R_plasma, and V_plasma. These are all
% vectors with the same length as the number of entries in the time vector.
% OUTPUTS
% energies - a structure holding the energies in the various components
% over the simulation run. This is in Mega-joules

mu0=4*pi*1e-7;
eV=11605;
Ru=8314.5;
Rgas=Ru/plasmaInps.MW; %this requires temeprature in KELVIN

% ---------- Energies in the circuit
MJ = 1e-6;
% E_fcc =  .5 *mu0 * Nfcc^2  * dfun(R_plasma)/2 .* I1.^2 *MJ;
energies.E_fcc =  .5 *mu0 * circInps.Nfcc^2  * circInps.Rfcc/2 .* modelOutputs.I1.^2 *MJ;
energies.E_plasma = .5*plasmaInps.m_propellant*modelOutputs.V_plasma.^2 *MJ;
energies.E_cap = .5*circInps.C*modelOutputs.Vcap.^2 *MJ;

energies.E_R1=ones(length(modelOutputs.t),1);
for i=1:length(energies.E_R1)
    if i == 1
        energies.E_R1(i)=0;
    else
        energies.E_R1(i)=circInps.R1*trapz(modelOutputs.t(1:i),modelOutputs.I1(1:i).^2)*MJ;
    end
end

energies.E_R2=ones(length(modelOutputs.t),1);
for i=1:length(energies.E_R2)
    if i == 1
        energies.E_R2(1)=0;
    else
        energies.E_R2(i)=circInps.R2*trapz(modelOutputs.t(1:i),modelOutputs.I2(1:i).^2)*MJ;
    end
end
energies.E_ohmic_losses=energies.E_R1+energies.E_R2;

gamma=plasmaInps.g;

T_0=plasmaInps.T0*1e3*eV+0.5*(gamma-1)*(modelOutputs.V_plasma(1).^2/(gamma*Rgas));

Tf=T_0-0.5*(gamma-1)*(modelOutputs.V_plasma.^2/(gamma*Rgas));

energies.E_therm=3*plasmaInps.m_propellant*Rgas*Tf*MJ;

energies.E_tot_2 = energies.E_fcc + energies.E_plasma + energies.E_cap; %needed for graph display script

energies.E_tot = energies.E_tot_2 + energies.E_R1 + energies.E_R2 + energies.E_therm; %needed for graph display script
end
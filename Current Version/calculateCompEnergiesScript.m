%% Energy Calculation Script
% Nathan Schilling
% 05/05/2020
% Holds all the necessary ancilary calcualtions needed to calculate various
% energies in the circuit (Fcc energy, plasma KE, load capacitor PE,
% resisitve losses, thermal energy using isentropic flow assumption)


% ---------- Energies in the circuit
MJ = 1e-6;
% E_fcc =  .5 *mu0 * Nfcc^2  * dfun(R_plasma)/2 .* I1.^2 *MJ;
E_fcc =  .5 *mu0 * Nfcc^2  * Rfcc/2 .* I1.^2 *MJ;
E_plasma = .5*m_propellant*V_plasma.^2 *MJ;
E_cap = .5*C*Vcap.^2 *MJ;

E_R1=ones(length(t),1);
for i=1:length(energies.E_R1)
    if i == 1
        E_R1(i)=0;
    else
        E_R1(i)=R1*trapz(t(1:i),I1(1:i).^2)*MJ;
    end
end

E_R2=ones(length(t),1);
for i=1:length(E_R2)
    if i == 1
        E_R2(1)=0;
    else
        E_R2(i)=R2*trapz(t(1:i),I2(1:i).^2)*MJ;
    end
end
E_ohmic_losses=E_R1+E_R2;

T_0=T0+0.5*(gamma-1)*(V_plasma(1).^2/(gamma*Rgas));

Tf=T_0-0.5*(gamma-1)*(V_plasma.^2/(gamma*Rgas));

energies.E_therm=3*m_propellant*Rgas*Tf*MJ;

energies.E_tot = E_fcc + E_plasma + E_cap + E_R1 + E_R2 + E_therm; %needed for graph display script

energies.E_tot_2 = E_fcc + E_plasma + E_cap; %needed for graph display script
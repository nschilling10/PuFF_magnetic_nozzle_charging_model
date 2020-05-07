function Pbrexact = calcRad(circInps,plasmaInps,modelOutputs)
%UNTITLED6 Calculates power radiated per sterradian per unit time over time 
%   Calculates power radiated per sterradian per unit time over time given
%   circuit parameters, parameters regarding the plasma, and the outputs 
%   from a simulation run

    %% Constants
    mu0=4*pi*1e-7;
    eV=11605;
    Ru=8314.5;
    %% Load in Plasma Parameters
    MW=plasmaInps.MW;
    Z=plasmaInps.Z;
    gamma=plasmaInps.g;
    Vp0=modelOutputs.V_plasma(1);
    %% Perform Calculations
    Rgas=Ru/MW; %this requires temeprature in KELVIN
    vol=pi*modelOutputs.R_plasma.^2.*modelOutputs.R_plasma; %calculating volume assuming plasma is cylinder with length scale equal to it's radius
    rho=plasmaInps.m_propellant./vol;
    T_0=plasmaInps.T0*1e3*eV+0.5*(gamma-1)*(Vp0.^2/(gamma*Rgas));
    Tf=T_0-0.5*(gamma-1)*(modelOutputs.V_plasma.^2/(gamma*Rgas));
 %   Pbrexact=51.6642e12.*sqrt(Tf).*Z.^3.*rho./MW.^2/(4*pi);
    Pbrexact=51.6642e12.*sqrt(plasmaInps.T0*1e3*eV).*Z.^3.*rho./MW.^2/(4*pi);
end


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
PUFF=true;
reactorGradeCaps=false;

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
    Rgas = 8314.5/MW;
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
    Rgas = 8314.5/plasmaInps.MW;
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
    
    plasmaInps.m_propellant = .155;
    %primary circuit (1) inputs
    circInps.L1 = 100e-9;
    circInps.R1 = 1e-5 * 0;
    
    %secondary circuit (2) inputs including reactor capacitance C
    circInps.L2 = 10e-9;
    circInps.R2 = 5e-3 * 0;
    circInps.C = .01;  %because we need 50 MJ or more
    
    %plasma inputs
    plasmaInps.MW = 235;
    Rgas = 8314.5/plasmaInps.MW;
    plasmaInps.T0 = 1;
    plasmaInps.g = 1.3;
    
    %transformer inputs
    circInps.mu_r = 1;
    circInps.k=0.9;
    
    circInps.rT = .1;
    circInps.NT1 = 100;
    circInps.NT2 = 1; %NT1/10;
    circInps.lT1 = circInps.rT*1;
    circInps.lT2 = circInps.rT*5;
    
    %flux compression generator inputs
    circInps.Nfcc = 1;
    circInps.Rfcc = 10;
    
    %-----------initial conditions
    circInps.I10 = 1e6; %Amps, initial current in primary
    circInps.I20 = 0;
%     circInps.V10 = 0;
    circInps.V20 = 0;  %voltage on reactor caps
end
%%
graphDisplay=true;
chargingModelOutputs=charging_nozzle_model(circInps,plasmaInps,graphDisplay);
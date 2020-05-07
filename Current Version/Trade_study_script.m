%% Trade study script
% study trades NT1 and NT2 to maximize the voltage on the capacitors. This
% will set NT1 and NT2 for a given target V_plasma
% Nathan Schilling
% 05/05/2020
clear all
close all

NT1_vec=1:1e2;
NT2_vec=1:1e2;

plasmaInps.m_propellant = .125;
%primary circuit (1) inputs
circInps.L1 = 100e-9;
circInps.R1 = 1e-5 * 0;

%secondary circuit (2) inputs including reactor capacitance C
circInps.L2 = 10e-9;
circInps.R2 = 5e-3 * 0;
circInps.C = .01;  %because we need 50 MJ or more

%plasma inputs
plasmaInps.MW = 4;
Rgas = 8314.5/plasmaInps.MW;
plasmaInps.T0 = 1; %[T0]=keV
plasmaInps.g = 1.3;

%transformer inputs
circInps.mu_r = 1;
circInps.k=0.9;

circInps.rT = .1;
circInps.lT1 = circInps.rT*1;
circInps.lT2 = circInps.rT*1;
        
%flux compression generator inputs
circInps.Nfcc = 100;
circInps.Rfcc = 1.7;

%-----------initial conditions
circInps.I10 = 1e6; %Amps, initial current in primary
circInps.I20 = 0;
circInps.V10 = 0;
circInps.V20 = 0;  %voltage on reactor caps

%------Energy the capacitors must be charged to in order for the Z-machine to work
targetEnergy=10*1e6;

graphDisplay=false;

Vpf_mat=zeros(length(NT1_vec),length(NT2_vec));
Vcapf_mat=zeros(length(NT1_vec),length(NT2_vec));

for i=1:length(NT1_vec)
    for j=1:length(NT1_vec)
        
        circInps.NT1 = NT1_vec(i);
        circInps.NT2 = NT2_vec(j);
        
        modelOutput=charging_nozzle_model(circInps,plasmaInps,graphDisplay); %run the model
        
        Vpf_mat(i,j)=-modelOutput.V_plasma(end);
        Vcapf_mat(i,j)=modelOutput.Vcap(end);
    end
end
%% Visualization
figure(11)
surf(NT2_vec,NT1_vec,Vcapf_mat*1e-3)
set(gca,'xscale','log')
set(gca,'yscale','log')
set(gca,'zscale','log')
xlabel('$N_{T2}$ (\#)','interpreter','latex','fontsize',24)
ylabel('$N_{T1}$ (\#)','interpreter','latex','fontsize',24)
zlabel('\boldmath$V_{cap_{f}}$ (kV)','interpreter','latex','fontsize',24)
figure(10)
surf(NT2_vec,NT1_vec,Vpf_mat*1e-3)
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('$N_{T2}$ (\#)','interpreter','latex','fontsize',24)
ylabel('$N_{T1}$ (\#)','interpreter','latex','fontsize',24)
zlabel('\boldmath$V_{plasma_{f}}$ (km/s)','interpreter','latex','fontsize',24)
%% Find max
[cols,rowInd_vec]=max(Vcapf_mat,[],1);
[val,colInd]=max(cols);
rowInd=rowInd_vec(colInd);
disp('NT1 best =')
NT1_vec(rowInd)
disp('Nt2 best =')
NT2_vec(colInd)
val == Vcapf_mat(rowInd,colInd)
circInps.NT1 = NT1_vec(rowInd);
circInps.NT1 = NT2_vec(colInd);

graphDisplay=true;
modelOutput=charging_nozzle_model(circInps,plasmaInps,graphDisplay); %run the model
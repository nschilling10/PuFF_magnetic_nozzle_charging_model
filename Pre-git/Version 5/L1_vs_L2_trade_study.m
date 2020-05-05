%% Nfcc vs. Rfcc trade study study
% study trades number of turns in the magnetic nozzle and the size of the
% magnetic nozzle to maxiize the exit velocity, and thereby the thrust and
% specific impulse.
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
circInps.lT2 = circInps.rT*4;

circInps.NT1 = 2;
circInps.NT2 = 20;

%flux compression generator inputs
circInps.Nfcc = 10;
circInps.Rfcc = 3.4;

%-----------initial conditions
circInps.I10 = 1e6; %Amps, initial current in primary
circInps.I20 = 0;
circInps.V10 = 0;
circInps.V20 = 0;  %voltage on reactor caps

graphDisplay=false;

EtaKE_mat=zeros(length(NT1_vec),length(NT2_vec));
EtaCirc_mat=zeros(length(NT1_vec),length(NT2_vec));
EtaSys_mat=zeros(length(NT1_vec),length(NT2_vec));

for i=1:length(NT1_vec)
    for j=1:length(NT2_vec)
        circInps.NT1 = NT1_vec(i);
        circInps.NT2 = NT2_vec(j);
        
        [eta_KE,eta_circ,eta_sys]=rlc_recharge_nozzle_circuit_v3_5(circInps,plasmaInps,graphDisplay);
        
        EtaKE_mat(i,j)=eta_KE;
        EtaCirc_mat(i,j)=eta_circ;
        EtaSys_mat(i,j)=eta_sys;
    end
end
%% Visualization
figure(10)
surf(NT1_vec,NT2_vec,EtaSys_mat)
set(gca,'xscale','log')
set(gca,'yscale','log')
set(gca,'zscale','log')
xlabel('$N_{T1}$ (\#)','interpreter','latex','fontsize',24)
ylabel('$N_{T2}$ (\#)','interpreter','latex','fontsize',24)
zlabel('\boldmath$\eta_{sys}$','interpreter','latex','fontsize',24)
%% Find max
[rows,colInd_vec]=max(EtaSys_mat);
[val,rowInd]=max(rows);
colInd=colInd_vec(rowInd);

circInps.NT1 = NT1_vec(colInd);
circInps.NT2 = NT2_vec(rowInd);

graphDisplay=true;

[eta_KE,eta_circ,eta_sys]= rlc_recharge_nozzle_circuit_v3_5(circInps,plasmaInps,graphDisplay);
eta_sys == val
%% Graph for KE efficiency
figure(11)
surf(NT1_vec,NT2_vec,EtaKE_mat)
set(gca,'xscale','log')
set(gca,'yscale','log')
set(gca,'zscale','log')
xlabel('$N_{T1}$ (\#)','interpreter','latex','fontsize',24)
ylabel('$N_{T2}$ (\#)','interpreter','latex','fontsize',24)
zlabel('\boldmath$\eta_{KE}$','interpreter','latex','fontsize',24)
[rows,colInd_vec]=max(EtaSys_mat);
[val,rowInd]=max(rows);
colInd=colInd_vec(rowInd);
NT1_vec(colInd);
NT2_vec(rowInd);
%% Graph for Circuit efficiency
figure(12)
surf(NT1_vec,NT2_vec,EtaCirc_mat)
set(gca,'xscale','log')
set(gca,'yscale','log')
set(gca,'zscale','log')
xlabel('$N_{T1}$ (\#)','interpreter','latex','fontsize',24)
ylabel('$N_{T2}$ (\#)','interpreter','latex','fontsize',24)
zlabel('\boldmath$\eta_{circ}$','interpreter','latex','fontsize',24)
[rows,colInd_vec]=max(EtaSys_mat);
[val,rowInd]=max(rows);
colInd=colInd_vec(rowInd);
NT1_vec(colInd);
NT2_vec(rowInd);
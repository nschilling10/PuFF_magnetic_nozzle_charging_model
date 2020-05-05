%% Nfcc vs. Rfcc trade study study
% study trades number of turns in the magnetic nozzle and the size of the
% magnetic nozzle to maxiize the exit velocity, and thereby the thrust and
% specific impulse.
% Nathan Schilling
% 05/05/2020
clear all
close all

Nfcc_vec=1:1e2;
Rfcc_vec=1:0.1:10;

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

circInps.NT1 = 2;
circInps.NT2 = 20;

%-----------initial conditions
circInps.I10 = 1e6; %Amps, initial current in primary
circInps.I20 = 0;
circInps.V10 = 0;
circInps.V20 = 0;  %voltage on reactor caps

graphDisplay=false;

Vpf_mat=zeros(length(Nfcc_vec),length(Rfcc_vec));

for i=1:length(Nfcc_vec)
    for j=1:length(Rfcc_vec)
        
        %flux compression generator inputs
        circInps.Nfcc = Nfcc_vec(i);
        circInps.Rfcc = Rfcc_vec(j);
        
        [Vp_f,capIsCharged]=rlc_recharge_nozzle_circuit_v3_8_1(circInps,plasmaInps,graphDisplay);
        
        if capIsCharged
            Vpf_mat(i,j)=Vp_f;
        else
            Vpf_mat(i,j)=NaN;
        end
    end
end
%% Data conditioning
[Nfcc_mat,Rfcc_mat]=meshgrid(Nfcc_vec,Rfcc_vec);
compMat=isfinite(Vpf_mat);
Nfcc_vec=Nfcc_mat(compMat);
Rfcc_vec=Rfcc_mat(compMat);
Vpf_vec=Vpf_mat(compMat);
%% Visualization
figure(10)
plot3(Nfcc_vec,Rfcc_vec,Vpf_vec*1e-3)
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('$N_{fcc}$ (\#)','interpreter','latex','fontsize',24)
ylabel('$R_{fcc}$ (\#)','interpreter','latex','fontsize',24)
zlabel('\boldmath$V_{plasma_{f}}$ (km/s)','interpreter','latex','fontsize',24)
[rows,colInd_vec]=max(EtaSys_mat);
[val,rowInd]=max(rows);
colInd=colInd_vec(rowInd);
NT1_vec(colInd);
NT2_vec(rowInd);
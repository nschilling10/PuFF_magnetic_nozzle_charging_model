%% Primary vs. secondary transformer winding study
% study trades primary and secondary number of turns on the transformer
% to find the optimal ratio for minimized current in the seed coils,
% provided the capacitor reaches it's charging voltage.
% Nathan Schilling
% 05/01/2020
clear all
close all

MA=1e6;

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

%flux compression generator inputs
circInps.Nfcc = 10;
circInps.Rfcc = 3.4;

%-----------initial conditions
circInps.I10 = 1e6; %Amps, initial current in primary
circInps.I20 = 0;
circInps.V10 = 0;
circInps.V20 = 0;  %voltage on reactor caps

graphDisplay=false;

Min_current_ratio_mat=ones(length(NT1_vec),length(NT2_vec));

for i=1:length(NT1_vec)
    for j=1:length(NT2_vec)
        circInps.NT1 = NT1_vec(i);
        circInps.NT2 = NT2_vec(j);
        
        [I_max,capIsCharged]=rlc_recharge_nozzle_circuit_v3_7(circInps,plasmaInps,graphDisplay);
        
        if capIsCharged
            Min_current_ratio_mat(i,j)=I_max/circInps.I10;
        end
    end
end
%% Visualization
figure(10)
surf(NT1_vec,NT2_vec,Min_current_ratio_mat)
set(gca,'xscale','log')
set(gca,'yscale','log')
set(gca,'zscale','log')
xlabel('$N_{T1}$ (\#)','interpreter','latex','fontsize',24)
ylabel('$N_{T2}$ (\#)','interpreter','latex','fontsize',24)
zlabel('\boldmath$\frac{I_{1_{min}}}{I_{1_{0}}}$','interpreter','latex','fontsize',24)
%% Find min
[rows,colInd_vec]=min(Min_current_ratio_mat);
[val,rowInd]=min(rows);
colInd=colInd_vec(rowInd);

circInps.NT1 = NT1_vec(rowInd);
circInps.NT2 = NT2_vec(colInd);

graphDisplay=true;

[I_max,capIsCharged]=rlc_recharge_nozzle_circuit_v3_7(circInps,plasmaInps,graphDisplay);
I_max/circInps.I10 == val
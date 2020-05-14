%% Trade study script
% study trades Nfcc, Rfcc, and I10 to maximize the final plasma velocity. This
% will set Nfcc, Rfcc, and I10 for a given target V_plasma
% Nathan Schilling
% 05/11/2020
clear all
close all

Nfcc_vec=round(logspace(1,3,50),0);
Rfcc_vec=round(logspace(-1,2,50),2);
I10_vec=logspace(-1,1,50)*1e6;

plasmaInps.m_propellant = 150e-3;
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
plasmaInps.T0 = 1; %[T0]=keV
plasmaInps.g = 1.3;

%transformer inputs
circInps.mu_r = 1;
circInps.k=0.9;

circInps.rT = .1;
circInps.lT1 = circInps.rT*1;
circInps.lT2 = circInps.rT*1;

circInps.NT1 = 300;
circInps.NT2 = 1;


%-----------initial conditions
circInps.I20 = 0;
circInps.V10 = 0;
circInps.V20 = 0;  %voltage on reactor caps
plasmaInps.Rp0 = 1e-2; 

%------Energy the capacitors must be charged to in order for the Z-machine to work
targetEnergy=50*1e6;

graphDisplay=false;

Vpf_tensor=zeros(length(Nfcc_vec),length(Rfcc_vec),length(I10_vec));
Vcapf_tensor=zeros(length(Nfcc_vec),length(Rfcc_vec),length(I10_vec));

for i=1:length(Nfcc_vec)
    for j=1:length(Rfcc_vec)
        for k=1:length(I10_vec)
            
            %flux compression generator inputs
            circInps.Nfcc = Nfcc_vec(i);
            circInps.Rfcc = Rfcc_vec(j);
            
            circInps.I10 = I10_vec(k); %Amps, initial current in primary
            
            modelOutput=charging_nozzle_model(circInps,plasmaInps,graphDisplay); %run the model
            
            Vpf_tensor(i,j,k)=-modelOutput.V_plasma(end);
            Vcapf_tensor(i,j,k)=modelOutput.Vcap(end);
        end
    end
end
%% Find max
[max_mat,page_mat]=max(Vcapf_tensor,[],3);
[cols,rowInd_vec]=max(max_mat,[],1);
[val,colInd]=max(cols);
rowInd=rowInd_vec(colInd);
pageInd=page_mat(rowInd,colInd);
disp('Nfcc best =')
Nfcc_vec(rowInd)
disp('Rfcc best =')
Rfcc_vec(colInd)
disp('I10 best =')
I10_vec(pageInd)
if val ~= Vcapf_tensor(rowInd,colInd,pageInd)
    disp('Indicies do not match max value')
    keyboard
end
circInps.Nfcc = Nfcc_vec(rowInd);
circInps.Rfcc = Rfcc_vec(colInd);
circInps.I10 = I10_vec(pageInd);

graphDisplay=true;
modelOutput=charging_nozzle_model(circInps,plasmaInps,graphDisplay); %run the model
%% Visualization
figure(11)
surf(Rfcc_vec,Nfcc_vec,Vcapf_tensor(:,:,pageInd)*1e-3)
set(gca,'xscale','log')
set(gca,'yscale','log')
set(gca,'zscale','log')
xlabel('$R_{fcc}$ (\#)','interpreter','latex','fontsize',24)
ylabel('$N_{fcc}$ (\#)','interpreter','latex','fontsize',24)
zlabel('\boldmath$V_{cap_{f}}$ (kV)','interpreter','latex','fontsize',24)
title('$V_{cap_{f}}$ for $I_{1_{0_{max}}}$','interpreter','latex','fontsize',28)
figure(10)
surf(Rfcc_vec,Nfcc_vec,Vpf_tensor(:,:,pageInd)*1e-3)
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('$R_{fcc}$ (\#)','interpreter','latex','fontsize',24)
ylabel('$N_{fcc}$ (\#)','interpreter','latex','fontsize',24)
zlabel('\boldmath$V_{plasma_{f}}$ (km/s)','interpreter','latex','fontsize',24)
title('$V_{plasma_{f}}$ for $I_{1_{0_{max}}}$','interpreter','latex','fontsize',28)
%% Changing page index visualization
pageInd=1;
surf(Rfcc_vec,Nfcc_vec,Vcapf_tensor(:,:,pageInd)*1e-3)
set(gca,'xscale','log')
set(gca,'yscale','log')
set(gca,'zscale','log')
xlabel('$R_{fcc}$ (\#)','interpreter','latex','fontsize',24)
ylabel('$N_{fcc}$ (\#)','interpreter','latex','fontsize',24)
zlabel('\boldmath$V_{cap_{f}}$ (kV)','interpreter','latex','fontsize',24)
title(strcat('$V_{cap_{f}}$ for $I_{1_{0}}=$',num2str(I10_vec(pageInd))),'interpreter','latex','fontsize',28)
figure(10)
surf(Rfcc_vec,Nfcc_vec,Vpf_tensor(:,:,pageInd)*1e-3)
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('$R_{fcc}$ (\#)','interpreter','latex','fontsize',24)
ylabel('$N_{fcc}$ (\#)','interpreter','latex','fontsize',24)
zlabel('\boldmath$V_{plasma_{f}}$ (km/s)','interpreter','latex','fontsize',24)
title(strcat('$V_{plasma_{f}}$ for $I_{1_{0}}=$',num2str(I10_vec(pageInd))),'interpreter','latex','fontsize',28)
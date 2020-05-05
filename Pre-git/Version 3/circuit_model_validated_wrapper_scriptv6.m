%% Nathan Schilling
% Line resistances trade study
% 06/28/19
clear all
close all

test.L1=110e-6;
test.L2=2.2e-6;
test.L0=330e-6;
test.I0=1e6;
test.graphDisplay=false;

R1_vec=linspace(0,1,25);
R2_vec=linspace(0,1,25);
Gain_mat=ones(length(R1_vec),length(R2_vec));
Req_circ_e=ones(length(R1_vec),length(R2_vec));
for i=1:length(R2_vec)
    test.R2=R2_vec(i);
    for j=1:length(R1_vec)
        test.R1=R1_vec(j);
        [E_gain,E_circ] = circuitModelFunValidatedv2_0(test);
        Req_circ_e(j,i)=E_circ;
        Gain_mat(j,i)=E_gain/E_circ;
    end
end
%% Plot results
fontSize=24;
h=figure(11);
surf(R1_vec,R2_vec,Gain_mat)
set(gca,'zscale','log')
xlabel('\boldmath$R_1$\textup{, }\boldmath$\Omega$','interpreter','latex','fontsize',fontSize)
ylabel('\boldmath$R_2$\textup{, }\boldmath$\Omega$','interpreter','latex','fontsize',fontSize)
zlabel('\textbf{Gain}','interpreter','latex','fontsize',fontSize)
h.Children.LineWidth=2;
h.Children.FontSize=18;
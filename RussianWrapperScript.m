%% Nathan Schilling
% Russian Transformer paper circuit model wrapper script
clear all
close all
hold off
% hold on

test.L1=25e-9;
test.L2=(6.8e-6);
test.I0=0.5e6;
test.L0=1.4e-6;
test.C_Load=5e-6;
% 
% test.graphDisplay=true;
% [E_gain,E_circ] = circuitModelFunv2_5(test)
% Gain=E_gain/E_circ
%[E_gain,E_circ] = circuitModelFunRussian(test)
%[E_gain,E_circ] = circuitModelFunv4(test)

%%
R=0;
l_0=10e-9;

input.l_1=l_0;
input.l_2=3e-6;
input.R1=R;
input.R2=R;
input.C=test.C_Load;

input.L1=test.L1;

input.L2=test.L2;
input.k=0.96;

input.tau=120e-6;
tSpan=linspace(0,input.tau,1e3);
dt=tSpan(2)-tSpan(1);
T=200e6;
input.v_exp_hand=@(t) dt;
input.dL_nozz_hand=@(t) test.L0*(-1/T);
input.L_nozz_hand=@(d) test.L0*(1-d/T);
input.I0=test.I0;
input.Rp0=0;
input.R_Fcc=0;
[t,I_1,I_2,V_Cap] = circuitModelFunction_V_Cassibry(input);
% Visualization
h=figure(1);
plot(t*10^6,I_1*1e-6,t*10^6,I_2*1e-6)
grid on
xlabel('\textbf{Time since ignition, }\boldmath$\mu$\textbf{s}','interpreter','latex','fontsize',22)
ylabel('\textbf{Current, MA}','interpreter','latex','fontsize',22)
title('\textbf{Current vs. Time since ignition with load connected}','interpreter','latex','fontsize',22)
legend({'\boldmath $I_1$','\boldmath $I_2$'},'interpreter','latex','fontsize',18)
h.Children(2).LineWidth=2;
h.Children(2).FontSize=18;
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
xlim([0 10])

figure(2)
plot(t*10^6,V_Cap*1e-3)
grid on
title('\textbf{U(kV) vs. Time since ignition}','interpreter','latex')
xlabel('\textbf{Time since ignition (}{\boldmath$\mu$}\textbf{sec)}','interpreter','latex')
ylabel('\textbf{Voltage across capacitor \textup{U}(kV)}','interpreter','latex')
set(gca,'fontsize',28)
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
xlim([0 10])

%%

figure(3)
legend('5\muC','50\muC','500\muC')

figure(4)
legend('500\muC','50\muC','5\muC')
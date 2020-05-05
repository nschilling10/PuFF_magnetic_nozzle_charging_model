%% Wrapper script
% Nathan Schilling
% 03/31/20
clear all
close all
% Universal constants
g0=9.8066; %[g0]=m/s^2
Ru=8314.3; %[Ru]=J/k-mol-K
k_b=8.617333262145e-5; %[k_b]=eV/K
e=1.602176634e-19; %[e]=J/eV
eps0=8.8541878128e-12; %[eps0]=SI units
mu_0=4*pi*1e-7; %[mu_0]=SI units
c=299792458; %[c]=m/s
h=6.62607015e-34; %[h]=J*s
m_a=1.6605e-27; %[m_a]=kg/u
m_e=9.10938356e-31; %[m_e]=kg
N_A=6.02214076e23; %[N_A]=#/mol
% Circuit input parameters
l_0=100e-9;
R=0;
C=4e-3;
input.l_1=l_0;
input.l_2=l_0;
input.R1=R;
input.R2=R;
input.C=C;
% Transformer params
mu_r=1;
r_T=0.1; %[r]=m
N_T1=25;
l_T1=1; %[l_T1]=m
input.L1=mu_0*mu_r*(pi*r_T^2)*N_T1.^2/l_T1;

N_T2=4;
l_T2=3; %[l_T1]=m
input.k=0.9;
input.L2=mu_0*mu_r*(pi*r_T^2)*N_T2.^2/l_T2;

% Seed coil parameters
N_Fcc=10;
r_Fcc=1; %[r_Fcc]=m
input.I0=1e4;

input.R_Fcc=r_Fcc;
L_Fcc=mu_0*(pi*r_Fcc^2)*N_Fcc^2/r_Fcc;
B_0=mu_0*mu_r*N_Fcc*input.I0/r_Fcc;
A_0=pi*r_Fcc^2;
Phi_0=L_Fcc*input.I0;
B=@(r_p) Phi_0./abs(A_0-pi*r_p.^2);
input.P_mag=@(r_p) B(r_p).^2/(2*mu_0);
%B=@(I1) mu_0*mu_r*N_Fcc*I1/r_Fcc;
%input.P_mag=@(I1) B(I1).^2/(2*mu_0);
input.V_sqig=@(r_p,v_p) 2*pi*B_0*r_p.*v_p./(1-(r_p/input.R_Fcc).^2);

% Plasma parameters
yield=10e6; %[yield]=J
input.m_p=0.05; %[m_p]=kg
Rp0=1e-2;
g=1.3;
MW=4; % need this

V_0=pi*Rp0^2*r_Fcc; % note this stands for the inital volume of the plasma -> plasma is cylindrical
rho_0=input.m_p/V_0;
n_0=rho_0*(N_A/MW);

T=yield/(e*k_b*n_0); %[T]=K

R_gas=Ru/MW;
%input.vp0=sqrt(g*R_gas*2*T);
input.vp0=sqrt(2*R_gas*T);
input.tau=(r_Fcc/input.vp0);
input.Rp0=0;

input.tauPerc=2;

[t,I_1,I_2,V_Cap,Rp_vec,Vp_vec] = circuitModelFunction_V_Cassibry2(input);
%% Post-processing

V_sqig_vec=input.V_sqig(Rp_vec,Vp_vec);

[V_Cap_maxVal,V_Cap_maxInd]=max(V_Cap);
E_cap=0.5*input.C*(V_Cap_maxVal^2-V_Cap(1)^2);

% Calulcate current changes
M_circ=input.k*sqrt(input.L1*input.L2);
dI1_vec=-(M_circ*V_Cap - input.L2*V_sqig_vec - input.l_2*V_sqig_vec ...
    + I_1*input.l_2*input.R1 + I_1*input.L2*input.R1 + I_2*M_circ*input.R2)...
            /(- M_circ^2 + input.l_1*input.l_2 + input.l_1*input.L2 + ...
            input.l_2*input.L1 + input.L1*input.L2);
dI2_vec=-(input.l_1*V_Cap + input.L1*V_Cap - M_circ*V_sqig_vec + I_2*input.l_1*input.R2...
        + I_2*input.L1*input.R2 + I_1*M_circ*input.R1)...
        /(- M_circ^2 + input.l_1*input.l_2 + input.l_1*input.L2 + ...
            input.l_2*input.L1 + input.L1*input.L2);
% Calculate energy of each component
E_Cap=0.5*input.C*V_Cap(1:V_Cap_maxInd).^2;

E_in=0.5*input.m_p*Vp_vec(V_Cap_maxInd)^2-yield;

E_gen=zeros(V_Cap_maxInd,1);
for i=1:length(E_gen)
    if i == 1
        E_gen(i)=0;
    else
        E_gen(i)=trapz(t(1:i),I_1(1:i).*V_sqig_vec(1:i));
    end
end

E_L1=0.5*input.L1*I_1(1:V_Cap_maxInd).^2;
E_L2=0.5*input.L2*I_2(1:V_Cap_maxInd).^2;
E_M=0.5*M_circ*(I_1(1:V_Cap_maxInd)-I_2(1:V_Cap_maxInd)).^2;

E_l1=0.5*input.l_1*I_1(1:V_Cap_maxInd).^2;
E_l2=0.5*input.l_2*I_2(1:V_Cap_maxInd).^2;

E_R1=ones(V_Cap_maxInd,1);
for i=1:length(E_R1)
    if i == 1
        E_R1(i)=0;
    else
        E_R1(i)=-trapz(t(1:i),input.R1*(I_1(1:i)).^2);
    end
end
E_R2=ones(V_Cap_maxInd,1);
for i=1:length(E_R2)
    if i == 1
        E_R2(i)=0;
    else
        E_R2(i)=-trapz(t(1:i),input.R2*(I_2(1:i)).^2);
    end
end

E_tot=E_Cap+E_in+E_gen+E_L1+E_L2+E_M+E_l1+E_l2+E_R1+E_R2;
Gain=E_cap/(E_tot(end)-E_tot(1)-E_R1(end)-E_R2(end))

Gain=E_cap/-E_in(end)
%% Plotting
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
%xlim([0 10])

figure(2)
plot(t*10^6,V_Cap*1e-3)
grid on
title('\textbf{U(kV) vs. Time since ignition}','interpreter','latex')
xlabel('\textbf{Time since ignition (}{\boldmath$\mu$}\textbf{sec)}','interpreter','latex')
ylabel('\textbf{Voltage across capacitor \textup{U}(kV)}','interpreter','latex')
set(gca,'fontsize',28)
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
%xlim([0 10])

% Domain of time from t=0 to t=V_cap_max
t_max_domain=t(1:V_Cap_maxInd);

figure(3)
plot(t_max_domain*10^6,E_tot*10^-6)
xlabel('Time since ignition (\mus)')
ylabel('Total energy in the circuit (MJ)')

figure(4)
plot(t_max_domain*10^6,E_Cap*10^-6,t_max_domain*10^6,E_in*10^-6,t_max_domain*10^6,E_gen*10^-6,t_max_domain*10^6,E_L1*10^-6,t_max_domain*10^6,E_L2*10^-6,t_max_domain*10^6,E_M*10^-6,t_max_domain*10^6,E_l1*10^-6,t_max_domain*10^6,E_l2*10^-6,t_max_domain*10^6,E_R1*10^-6,'r--',t_max_domain*10^6,E_R2*10^-6,'g--')
xlabel('Time since ignition (\mus)')
ylabel('Energy in component (MJ)')
legend('E_{Cap}','E_{in}','Generator','L_1','L_2', 'M','l_1','l_2','R_1','R_2','Location','Northeast')

figure(5)
plot(t_max_domain*10^6,E_Cap*10^-6,t_max_domain*10^6,E_in*10^-6,t_max_domain*10^6,E_gen*10^-6,t_max_domain*10^6,E_L1*10^-6,t_max_domain*10^6,E_L2*10^-6,t_max_domain*10^6,E_M*10^-6)
xlabel('Time since ignition (\mus)')
ylabel('Energy in component (MJ)')
legend('E_{Cap}','E_{in}','Generator','L_1','L_2', 'M','Location','Northeast')

figure(6)
plot(t*1e6,Rp_vec)
xlabel('Time since ignition (\mus)','fontsize',18)
ylabel('$R_{p}$ (m)','interpreter','latex','fontsize',20)

figure(7)
plot(t*1e6,Vp_vec*1e-3)
xlabel('Time since ignition (\mus)','fontsize',18)
ylabel('$v_{exp}$ (km/s)','interpreter','latex','fontsize',20)
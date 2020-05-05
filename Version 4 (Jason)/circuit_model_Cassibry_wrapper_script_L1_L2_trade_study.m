%% Problem 1b-c script
% Nathan Schilling
% 11/24/19
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
% Problem input parameters
l_0=100e-9;
R=1;
C=400e-6;
input.l_1=l_0;
input.l_2=l_0;
input.R1=R;
input.R2=R;
input.C=C;
% Gas params
T=1e3/k_b; %[T_i]=eV
g=1.3;
MW=4; % need this

% Transformer params
mu_r=1;
r_T=0.1; %[r]=m
N_T1_vec=logspace(0,3,1e2);
l_T1=1; %[l_T1]=m
L1_vec=mu_0*mu_r*(pi*r_T^2)*N_T1_vec.^2/l_T1;

N_T2_vec=logspace(0,3,1e2);
l_T2=3; %[l_T1]=m
L2_vec=mu_0*mu_r*(pi*r_T^2)*N_T2_vec.^2/l_T2;
input.k=0.9;

N_Fcc=10;
r_Fcc=3.4; %[r_Fcc]=m
L_Fcc=mu_0*(pi*r_Fcc^2)*N_Fcc^2/r_Fcc;
R_gas=Ru/MW;
input.tau=(r_Fcc/(4*sqrt(g*R_gas*T)));
input.v_exp_hand= @(t) sqrt(g*R_gas*2*T)*(cos(pi*t/(2*input.tau)).*...
    ((t/input.tau)<1)-(1.*(t/input.tau)>=Inf));
input.dL_nozz_hand=@(t) -0.5*mu_0*N_Fcc*input.v_exp_hand(t);
input.L_nozz_hand=@(d) 0.5*mu_0*N_Fcc*d;
input.I0=5e6;
input.Rp0=0;
input.R_Fcc=r_Fcc;

Gain_mat=ones(length(L2_vec),length(L1_vec));
for i=1:length(L1_vec)
    for j=1:length(L2_vec)
        input.L1=L1_vec(i);
        input.L2=L2_vec(j);
        [t,I_1,I_2,V_Cap,d_vec] = circuitModelFunction_V_Cassibry2(input);
        L_t_vec=input.L_nozz_hand(d_vec);
        [V_Cap_maxVal,V_Cap_maxInd]=max(V_Cap);
        E_gain=0.5*input.C*(V_Cap(V_Cap_maxInd)^2-V_Cap(1)^2);
        V_sqig=I_1.*input.dL_nozz_hand(t);
%     % Calulcate current changes
%     L0_vec=input.L_T1+input.L_1+L_t_vec;
%     V_sqig_vec=I_1.*input.dLt_handle(t);
%     dI2_dt_vec=(V_Cap-I_2*input.R2-(M_circ*(V_sqig_vec-I_1*input.R1)./L0_vec))...
%         ./(input.L_2+input.L_T2-(M_circ^2./L0_vec));
%     dI1_dt_vec=(V_sqig_vec-I_1*input.R1-M_circ*dI2_dt_vec)./(L0_vec);
%        E_circ=-trapz(t(1:V_Cap_maxInd),(I_1(1:V_Cap_maxInd)-input.I0).*V_sqig(1:V_Cap_maxInd));
        E_circ=-trapz(t(1:V_Cap_maxInd),I_1(1:V_Cap_maxInd).*V_sqig(1:V_Cap_maxInd));
        Gain_mat(j,i) = E_gain/E_circ;
    end
end
%% Plotting 
figure(11);
surf(N_T1_vec,N_T2_vec,Gain_mat)
set(gca,'xscale','log')
set(gca,'yscale','log')
set(gca,'zscale','log')
xlabel('\boldmath $N_{T1}$','interpreter','latex','fontsize',24)
ylabel('\boldmath $N_{T2}$','interpreter','latex','fontsize',24)
zlabel('\textbf{Gain (}{\boldmath$\frac{\Delta E_{cap}}{\Delta E_{gen}}$}\textbf{)}','interpreter','latex','fontsize',36)
%% Ideal Case plot
[rows,colInd_vec]=max(Gain_mat);
[val,rowInd]=max(rows);
colInd=colInd_vec(rowInd);
input.L1=L1_vec(rowInd);
input.L2=L2_vec(colInd);

[t,I_1,I_2,V_Cap] = circuitModelFunction_V_Cassibry2(input);

h=figure(1);
plot(t*10^6,I_1*1e-6,t*10^6,I_2*1e-6)
grid on
xlabel('\textbf{Time since ignition, }\boldmath$\mu$\textbf{s}','interpreter','latex','fontsize',22)
ylabel('\textbf{Current, MA}','interpreter','latex','fontsize',22)
%title('\textbf{Current vs. Time since ignition with load connected}','interpreter','latex','fontsize',22)
legend({'\boldmath $I_1$','\boldmath $I_2$'},'interpreter','latex','fontsize',18)
h.Children(2).LineWidth=2;
h.Children(2).FontSize=18;
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
xlim([0 10])

figure(2)
plot(t*10^6,V_Cap*1e-3)
grid on
%title('\textbf{U(kV) vs. Time since ignition}','interpreter','latex')
xlabel('\textbf{Time since ignition (}{\boldmath$\mu$}\textbf{sec)}','interpreter','latex')
ylabel('\textbf{Voltage across capacitor \textup{U}(kV)}','interpreter','latex')
set(gca,'fontsize',18)
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
xlim([0 10])

[V_Cap_maxVal,V_Cap_maxInd]=max(V_Cap);
E_gain=0.5*input.C*(V_Cap(V_Cap_maxInd)^2-V_Cap(1)^2);
V_sqig=I_1.*input.dL_nozz_hand(t);
%E_circ=-trapz(t(1:V_Cap_maxInd),(I_1(1:V_Cap_maxInd)-input.I0).*V_sqig(1:V_Cap_maxInd));
E_circ=-trapz(t,I_1.*V_sqig)
E_gain/E_circ
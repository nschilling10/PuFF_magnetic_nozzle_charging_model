%% Nathan Schilling
% Transformer circuit model wrapper script
% 02/19/19
clear all
close all
format long
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

% Circuit input parameters
l_0=100e-9;
R=10e-3;
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
N_T1=25;
l_T1=1; %[l_T1]=m
input.L1=mu_0*mu_r*(pi*r_T^2)*N_T1.^2/l_T1;

N_T2=4;
l_T2=3; %[l_T1]=m
input.L2=mu_0*mu_r*(pi*r_T^2)*N_T2.^2/l_T2;
input.k=0.9;

N_Fcc=10;
r_Fcc=3.4; %[r_Fcc]=m
L_Fcc=mu_0*(pi*r_Fcc^2)*N_Fcc^2/r_Fcc;
R_gas=Ru/MW;
input.tau=(r_Fcc/(4*sqrt(g*R_gas*T)));
input.v_exp_hand= @(t) sqrt(g*R_gas*2*T)*(cos(pi*t/(2*input.tau)).*...
    ((t/input.tau)<1)-(1.*(t/input.tau)>Inf));
input.dL_nozz_hand=@(t) -0.5*mu_0*N_Fcc*input.v_exp_hand(t);
input.L_nozz_hand=@(d) 0.5*mu_0*N_Fcc*d;
input.I0=5e6;
input.Rp0=0;
input.R_Fcc=r_Fcc;

[t_vec,I_1,I_2,V_Cap,d_vec] = circuitModelFunValidatedv2_1(input);

E_gain=0.5*input.C_load*(V_Cap(end)^2-V_Cap(1)^2);
E_in=0.5*trapz(t_vec,I_1.^2.*dL_nozz_hand(t_vec));
%% Visualization
visualization=true;
if visualization
    % Calulcate current changes
    L_01=input.L1+l_1+holderArray.Lt_handle(d);
    L_02=input.L2+l_2;
    M_circ=input.k*sqrt(input.L1*input.L2);
    
    dI1_dt=(M_circ*V_cap + I_1*L_02*input.R1 + I_2*M_circ*pinput.R2 + I_1c*L_02.*holderArray.dLt_handle(t_vec))./(M_circ^2 - L_01*L_02);
    dV_dt=I_2/input.C;
    dI2_dt=(L_01*V_cap + I_2*L_01*input.R2 + I_1*M_circ*R_1 + I_1*M_circ.*holderArray.dLt_handle(t_vec))./(M_circ^2 - L_01*L_02);
    dvec_dt=input.v_exp_hand(t_vec);
    
    figure(1)
    semilogy(t_vec*10^6,I_1,t_vec*10^6,I_2)
    grid on
    xlabel('\textbf{Time since ignition, }\boldmath$\mu$\textbf{s}','interpreter','latex','fontsize',22)
    ylabel('\textbf{Current, A}','interpreter','latex','fontsize',22)
    title('Current vs. Time since ignition with load connected')
    %legend({'\boldmath $I_1$','\boldmath $I_2$'},'interpreter','latex','fontsize',18)
    h.Children(2).LineWidth=2;
    h.Children(2).FontSize=18;
    set(findall(gca, 'Type', 'Line'),'LineWidth',2);
    
    %         fontSize=36;
    %         xlabel('\textbf{Time since ignition (}{\boldmath$\mu$}\textbf{sec)}','interpreter','latex','fontsize',fontSize)
    %         ylabel('\textbf{Current (A)}','interpreter','latex','fontsize',fontSize)
    
    h=figure(2);
    plot(t_vec*10^6,I_1*1e-6)
    set(findall(gca, 'Type', 'Line'),'LineWidth',2);
    grid on
    xlabel('\textbf{Time since ignition, }\boldmath$\mu$\textbf{s}','interpreter','latex','fontsize',22)
    ylabel('\boldmath$I_1$textbf{, MA}','interpreter','latex','fontsize',22)
    title('Generator Current vs. Time since ignition')
    h.Children.LineWidth=2;
    h.Children.FontSize=18;
    %         xlabel('\textbf{Time since ignition (}{\boldmath$\mu$}\textbf{sec)}','interpreter','latex','fontsize',fontSize)
    %         ylabel('\textbf{Current on generator side (MA)}','interpreter','latex','fontsize',fontSize)
    %         set(gca,'fontsize',28)
    %         set(gca,'LineWidth',2)
    
    figure(3)
    plot(t_vec*10^6,V_Cap*1e-3)
    grid on
    xlabel('Time since ignition (\mus)')
    ylabel('Capacitor voltage (kV)')
    title('Capacitor voltage vs. Time since ignition')
    %         xlabel('\textbf{Time since ignition (}{\boldmath$\mu$}\textbf{sec)}','interpreter','latex','fontsize',fontSize)
    %         ylabel('\textbf{Voltage across capacitor \textup{U}(kV)}','interpreter','latex','fontsize',fontSize)
    %         set(gca,'fontsize',28)
    %         set(gca,'LineWidth',2)
    
    figure(4)
    plot(t_vec*10^6,I_2*1e-6)
    grid on
    xlabel('Time since ignition (\mus)')
    ylabel('Load Current (MA)')
    title('Load Current vs. Time since ignition')
    %         xlabel('\textbf{Time since ignition (}{\boldmath$\mu$}\textbf{sec)}','interpreter','latex','fontsize',fontSize)
    %         ylabel('\textbf{Current on generator side (MA)}','interpreter','latex','fontsize',fontSize)
    %         set(gca,'fontsize',28)
    %         set(gca,'LineWidth',2)
    
    % Energy Calculation
    %        E_gen_ind=0.5*L_nozz_hand(t).*I_1.^2;
    E_L1=0.5*(input.L1-M_circ)*I_1.^2;
    E_l1=0.5*input.l_1*I_1.^2;
    E_l2=0.5*input.l_2*I_2.^2;
    E_L2=0.5*(input.L2-M_circ)*I_2.^2;
    E_M=0.5*M_circ*(I_1-I_2).^2;
    E_cap=0.5*input.C*V_Cap.^2;
    E_R1=ones(length(t_vec),1);
    for i=1:length(E_R1)
        if i == 1
            E_R1(i)=0;
        else
            E_R1(i)=-input.R1*trapz(t_vec(1:i),I_1(1:i).^2);
        end
    end
    E_R2=ones(length(t_vec),1);
    for i=1:length(E_R2)
        if i == 1
            E_R2(1)=0;
        else
            E_R2(i)=-input.R2*trapz(t_vec(1:i),I_2(1:i).^2);
        end
    end
    E_in=ones(length(t_vec),1);
    for i=1:length(E_in)
        if i == 1
            E_in(1)=0;
        else
            E_in(i)=trapz(t_vec(1:i),(I_1(1:i).^2.*dL_nozz_hand(t_vec(1:i))));
        end
    end
    E_coil=ones(length(t_vec),1);
    for i=1:length(E_coil)
        if i == 1
            E_coil(1)=0;
        else
            E_coil(i)=trapz(t_vec(1:i),(I_1(1:i).*L_nozz_hand(t_vec(1:i)).*dI1_dt(1:i)));
        end
    end
    E_gen=E_in+E_coil;
    E_tot=E_gen+E_L1+E_l1+E_l2+E_L2+E_M+E_cap+E_R1+E_R2;
    %         figure(5)
    %         plot(t,E_gen)
    %
    %         subplot(4,1,1), plot(t,E_gen)
    %         subplot(4,1,2), plot(t,E_L1)
    %         subplot(4,1,3), plot(t,E_L2)
    %         subplot(4,1,4), plot(t,E_M)
    %         figure(6)
    %         plot(t,E_cap)
    %         figure(7)
    %         subplot(2,1,1), plot(t,E_R1)
    %         subplot(2,1,2), plot(t,E_R2)
    %         figure(8)
    %         plot(t,E_tot)
    figure(5)
    plot(t_vec*10^6,E_in*10^-6,t_vec*10^6,E_coil*10^-6,t_vec*10^6,E_L1*10^-6,t_vec*10^6,E_L2*10^-6,t_vec*10^6,E_l1*10^-6,t_vec*10^6,E_l2*10^-6,t_vec*10^6,E_M*10^-6,t_vec*10^6,E_R1*10^-6,t_vec*10^6,E_R2*10^-6)
    xlabel('Time since ignition (\mus)')
    ylabel('Energy in component (MJ)')
    legend('E_{in}','E_{coil}','L_1','L_2','l_1','l_2','M','R_1','R_2','Location','Southwest')
    figure(6)
    plot(t_vec*10^6,(E_tot-E_tot(1))*10^-6)
    xlabel('Time since ignition (\mus)')
    ylabel('Total energy in the circuit (MJ)')
    
    % Voltage Calculation
    v_suigt=dL_nozz_hand(t_vec).*I_1;
    
    voltage_L1=input.L1.*dI1_dt;
    voltage_L2=input.L2.*dI2_dt;
    
    V_loop1= v_suigt+L_nozz_hand(t_vec).*dI1_dt+I_1*input.R1+voltage_L1 - M_circ*dI2_dt+input.l_1*dI1_dt;
    V_loop2= voltage_L2+I_2*input.R2+V_Cap - M_circ*dI1_dt+input.l_2*d21_dt;
    
    figure(7)
    plot(t_vec*1e6,V_loop1,t_vec*1e6,V_loop2)
    xlabel('Time since ignition (\mus)')
    ylabel('Voltage in a loop')
    legend('V loop 1','V loop 2','Generator voltage')
    
    figure(8)
    plot(t_vec*1e6,voltage_L1.*I_1,t_vec*1e6,voltage_L2.*I_2)
    xlabel('Time since ignition (\mus)')
    ylabel('Power')
    
    Power_vec_L1=abs(voltage_L1.*I_1);
    Power_vec_L2=abs(voltage_L2.*I_2);
    E_heat_L1=zeros(length(Power_vec_L1),1);
    E_heat_L2=zeros(length(Power_vec_L2),1);
    for i=2:length(Power_vec_L1)
        E_heat_L1(i)=trapz(t_vec(1:i),Power_vec_L1(1:i));
        E_heat_L2(i)=trapz(t_vec(1:i),Power_vec_L2(1:i));
    end
    max_pwr=max(voltage_L1.*I_1)
    max_pwr=max(voltage_L2.*I_2)
    max_heat_energy=max(max(E_heat_L1),max(E_heat_L2))
    figure(9)
    plot(t_vec*1e6,E_heat_L1,t_vec*1e6,E_heat_L2)
    xlabel('Time since ignition (\mus)')
    ylabel('Heating energy')
    figure(6)
        
end
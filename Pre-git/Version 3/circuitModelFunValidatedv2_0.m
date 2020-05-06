% Nathan Schilling
% 02/19/19
function [E_gain,E_circ] = circuitModelFunValidatedv2_0(circuitInputParams)
   
    % Physics constants
    mu0=4*pi*10^-7;
    
    % Nozzle Input parameters (SI units)
    x_0=1;
    A=pi*x_0^2;
    N=1e1;
    l_nozzle=1;
    v_0=(25e3/3);
    
    % Circuit parameters
    input.L1=circuitInputParams.L1;
    input.L2=circuitInputParams.L2;
    input.k=0.85;
    M_circ=input.k*sqrt(input.L1*input.L2);
    input.l_0=1e-9;
    if ~isfield(circuitInputParams,'R1')
        input.R1=0;
    else
        input.R1=circuitInputParams.R1;
    end
    if ~isfield(circuitInputParams,'R2')
        input.R2=0;
    else
        input.R2=circuitInputParams.R2;
    end
    input.C_load=22e-3;
    
    % Display parameters
    tauPerc=1;
    
    % Seed current
    I1_0=circuitInputParams.I0;
    
    % Nozzle inductive derived parameters
    global L0
    %L0=mu0*N^2*A/l_nozzle;
    L0=circuitInputParams.L0;
    global tau
    %tau=x_0/v_0;
    tau=10e-6;
    
    function out=L_nozz_hand(x)
        tau_1=2e-6;
%         if x < tau
            out=L0*(exp(-x.^2/(tau_1^2))+0.5*(tanh((x-tau_1*3)/tau_1*3)+1));
%         else
%             out=0;
%         end
    end

    function out=dL_nozz_hand(x)
        tau_1=2e-6;
%         if x < tau
            out=-L0*((3*(tanh((3*(3*tau_1 - x))/tau_1).^2 - 1))/(2*tau_1) + ...
                (2*x.*exp(-x.^2/tau_1^2))/tau_1^2);
%         else
%             out=0;
%         end
    end

    % Set up input functions
    D2L_nozz_hand=@(x) 0*x;
    input.Lt_handle=@(x) L_nozz_hand(x);
    input.dLt_handle=@(x) dL_nozz_hand(x);
    
    % Inital Condition
    Iode=[I1_0;0;0];
    
    tSpan=linspace(0,tauPerc*tau,1e3);
    
    %options=odeset('RelTol',1e-20,'AbsTol',1e-20);
    
    [t,y]=ode45(@(t,y) fcgfuns(t,y,input),tSpan,Iode);   

    %length(t)
    
    I_1=y(:,1);
    V_Cap=y(:,2);
    I_2=y(:,3);
    
    % Calulcate current changes
    L0_vec=input.L1+input.l_0+input.Lt_handle(t);
    dI1_dt=(M_circ*V_Cap + I_1*input.L2*input.R1 + I_2*M_circ*input.R2 + I_1.*input.L2.*input.dLt_handle(t))./(M_circ^2 - L0_vec*input.L2);
    dIV_dt=I_2/input.C_load;
    dI2_dt=(L0_vec.*V_Cap + I_2.*L0_vec*input.R2 + I_1*M_circ*input.R1 + I_1.*M_circ.*input.dLt_handle(t))./(M_circ^2 - L0_vec*input.L2);
    
    if circuitInputParams.graphDisplay
        
        figure(1)
        semilogy(t*10^6,I_1,t*10^6,I_2)
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
        plot(t*10^6,I_1*1e-6)
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
        plot(t*10^6,V_Cap*1e-3)
        grid on
        xlabel('Time since ignition (\mus)')
        ylabel('Capacitor voltage (kV)')
        title('Capacitor voltage vs. Time since ignition')
        %         xlabel('\textbf{Time since ignition (}{\boldmath$\mu$}\textbf{sec)}','interpreter','latex','fontsize',fontSize)
        %         ylabel('\textbf{Voltage across capacitor \textup{U}(kV)}','interpreter','latex','fontsize',fontSize)
        %         set(gca,'fontsize',28)
        %         set(gca,'LineWidth',2)
        
        figure(4)
        plot(t*10^6,I_2*1e-6)
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
        E_l=0.5*input.l_0*I_1.^2;
        E_L2=0.5*(input.L2-M_circ)*I_2.^2;
        E_M=0.5*M_circ*(I_1-I_2).^2;
        E_cap=0.5*input.C_load*V_Cap.^2;
        E_R1=ones(length(t),1);
        for i=1:length(E_R1)
            if i == 1
                E_R1(i)=0;
            else
                E_R1(i)=-input.R1*trapz(t(1:i),I_1(1:i).^2);
            end
        end
        E_R2=ones(length(t),1);
        for i=1:length(E_R2)
            if i == 1
                E_R2(1)=0;
            else
                E_R2(i)=-input.R2*trapz(t(1:i),I_2(1:i).^2);
            end
        end
%         E_gen=ones(length(t),1);
%         for i=1:length(E_gen)
%             if i == 1
%                 E_gen(1)=0;
%             else
%                 E_gen(i)=trapz(t(1:i),(I_1(1:i).*L_nozz_hand(t(1:i)).*dI1_dt(1:i)+I_1(1:i).^2.*dL_nozz_hand(t(1:i))));
%             end
%         end
        E_in=ones(length(t),1);
        for i=1:length(E_in)
            if i == 1
                E_in(1)=0;
            else
                E_in(i)=trapz(t(1:i),(I_1(1:i).^2.*dL_nozz_hand(t(1:i))));
            end
        end
        E_gen=ones(length(t),1);
        for i=1:length(E_gen)
            if i == 1
                E_gen(1)=0;
            else
                E_gen(i)=trapz(t(1:i),L_nozz_hand(t(1:i)).*I_1(1:i).*dI1_dt(1:i));
            end
        end
        E_tot=E_in+E_gen+E_L1+E_l+E_L2+E_M+E_cap+E_R1+E_R2;
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
        plot(t*10^6,(E_in-E_in(1))*10^-6,t*10^6,E_L1*10^-6,t*10^6,E_L2*10^-6,t*10^6,E_M*10^-6,t*10^6,E_R1*10^-6,t*10^6,E_R2*10^-6)
        xlabel('Time since ignition (\mus)')
        ylabel('Energy in component (MJ)')
        legend('Generator','L_1','L_2','M','R_1','R_2','Location','Southwest')
        figure(6)
        plot(t*10^6,(E_tot-E_tot(1))*10^-6)
        xlabel('Time since ignition (\mus)')
        ylabel('Total energy in the circuit (MJ)')
        
        % Voltage Calculation      
        voltage_in_circuit=dL_nozz_hand(t).*I_1;
        
        voltage_L1=input.L1.*dI1_dt;
        voltage_L2=input.L2.*dI2_dt;
        
        V_loop1= voltage_in_circuit+L_nozz_hand(t).*dI1_dt+I_1*input.R1+voltage_L1 - M_circ*dI2_dt+input.l_0*dI1_dt;
        V_loop2= voltage_L2+I_2*input.R2+V_Cap - M_circ*dI1_dt;
        
        figure(7)
        plot(t*1e6,V_loop1,t*1e6,V_loop2)
        xlabel('Time since ignition (\mus)')
        ylabel('Voltage in a loop')
        legend('V loop 1','V loop 2','Generator voltage')
        
        figure(8)
        plot(t*1e6,voltage_L1.*I_1,t*1e6,voltage_L2.*I_2)
        xlabel('Time since ignition (\mus)')
        ylabel('Power')
        
        Power_vec_L1=abs(voltage_L1.*I_1);
        Power_vec_L2=abs(voltage_L2.*I_2);
        E_heat_L1=zeros(length(Power_vec_L1),1);
        E_heat_L2=zeros(length(Power_vec_L2),1);
        for i=2:length(Power_vec_L1)
            E_heat_L1(i)=trapz(t(1:i),Power_vec_L1(1:i));
            E_heat_L2(i)=trapz(t(1:i),Power_vec_L2(1:i));
        end
        max_pwr=max(voltage_L1.*I_1)
        max_pwr=max(voltage_L2.*I_2)
        max_heat_energy=max(max(E_heat_L1),max(E_heat_L2))
        figure(9)
        plot(t*1e6,E_heat_L1,t*1e6,E_heat_L2)
        xlabel('Time since ignition (\mus)')
        ylabel('Heating energy')
        figure(6)
        
    end
    
    E_gain=0.5*input.C_load*(V_Cap(end)^2-V_Cap(1)^2);
    %E_circ=-trapz(t,I_1.*L_nozz_hand(t).*dI1_dt+I_1.^2.*dL_nozz_hand(t));
    E_circ=-trapz(t,I_1.^2.*dL_nozz_hand(t));
    
    %this is the ode circuit solver
    function dI = fcgfuns(t,y,holderArray)
        %varargin is a structure that is assumed to hold the following
        %fields
        I1=y(1);
        if ~isfield(holderArray,'Lt_handle')
            msgID = 'fcgfuns:BadInput';
            msg = 'No Lt_handle input field found.';
            baseException = MException(msgID,msg);
            throw(baseException);
        end
        
        if ~isfield(holderArray,'dLt_handle')
            msgID = 'fcgfuns:BadInput';
            msg = 'No dLt_handle input field found.';
            baseException = MException(msgID,msg);
            throw(baseException);
        end
        
        if ~isfield(holderArray,'L1')
            L_1=1;
        else
            L_1=holderArray.L1;
        end
        
        if ~isfield(holderArray,'L2')
            L_2=1;
        else
            L_2=holderArray.L2;
        end
        
        if ~isfield(holderArray,'k')
            M=0.85*sqrt(L_1*L_2);
        else
            M=holderArray.k*sqrt(L_1*L_2);
        end
        
        if ~isfield(holderArray,'l_0')
            l=0;
        else
            l=holderArray.l_0;
        end
        
        if ~isfield(holderArray,'R1')
            R_1=0;
        else
            R_1=holderArray.R1;
        end
        
        if ~isfield(holderArray,'R2')
            R_2=0;
        else
            R_2=holderArray.R2;
        end
        
        if ~isfield(holderArray,'C_load')
            C_load=1;
        else
            C_load=holderArray.C_load;
        end
        
        %L_0=holderArray.Lt_handle(t)+l+L_1-(L_1*L_2)/M;
        L_0=L_1+l+holderArray.Lt_handle(t);
        
        V_cap=y(2);
        I2=y(3);
        
        dI=zeros(3,1);
        dI(1,1)=(M*V_cap + I1*L_2*R_1 + I2*M*R_2 + I1*L_2*holderArray.dLt_handle(t))/(M^2 - L_0*L_2);
        dI(2,1)=I2/C_load;
        dI(3,1)=(L_0*V_cap + I2*L_0*R_2 + I1*M*R_1 + I1*M*holderArray.dLt_handle(t))/(M^2 - L_0*L_2);
    end
end
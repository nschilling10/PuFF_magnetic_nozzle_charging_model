% Nathan Schilling
% 02/19/19
function [E_gain,E_circ] = circuitModelFunv4(circuitInputParams)
   
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
    input.a=circuitInputParams.a;
    input.L2=(input.a)^2*input.L1;
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
    
    input.C_load=450e-3;
    
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
        if x < tau
            out=L0*(exp(-x.^2/(tau_1^2))+0.5*(tanh((x-tau_1*3)/tau_1*3)+1));
        else
            out=0;
        end
    end

    function out=dL_nozz_hand(x)
        tau_1=2e-6;
        if x < tau
            out=-L0*((3*(tanh((3*(3*tau_1 - x))/tau_1).^2 - 1))/(2*tau_1) + ...
                (2*x.*exp(-x.^2/tau_1^2))/tau_1^2);
        else
            out=0;
        end
    end

    % Set up input functions
    D2L_nozz_hand=@(x) 0*x;
    input.Lt_handle=@(x) L_nozz_hand(x);
    input.dLt_handle=@(x) dL_nozz_hand(x);
    
    % Inital Condition
    Iode=[I1_0;0;0];
    
    tSpan=linspace(0,tauPerc*tau,1e5);
    
    [t,y]=ode45(@(t,y) fcgfuns(t,y,input),tSpan,Iode);   

    I_1=y(:,1);
    V_Cap=y(:,2);
    I_2=y(:,3);
    
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

        E_gen=0.5*L_nozz_hand(t).*I_1.^2;
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
                E_R1(i)=input.R1*trapz(t(1:i),I_1(1:i).^2);
            end
        end
        E_R2=ones(length(t),1);
        for i=1:length(E_R2)
            if i == 1
                E_R2(1)=0;
            else
                E_R2(i)=input.R2*trapz(t(1:i),I_2(1:i).^2);
            end
        end
        E_tot=E_gen+E_L1+E_l+E_L2+E_M+E_cap-E_R1-E_R2;
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
    plot(t*10^6,E_gen*10^-6,t*10^6,E_L2*10^-6,t*10^6,E_M*10^-6,t*10^6,E_R1*10^-6,t*10^6,E_R2*10^-6)
    xlabel('Time since ignition (\mus)')
    ylabel('Energy in component (MJ)')
    legend('Generator','L_2','M','R_1','R_2','Location','Northwest')
    figure(6)
    plot(t*10^6,E_tot*10^-6)
    xlabel('Time since ignition (\mus)')
    ylabel('Total energy in the circuit (MJ)')

    figure(3)
    end
    
    E_gain=0.5*input.C_load*(V_Cap(end)^2-V_Cap(1)^2);
    E_circ=0.5*L0*I_1(1)^2;
    
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
        
        if ~isfield(holderArray,'a')
            a=1;
        else
            a=holderArray.a;
        end
        
        %L_0=holderArray.Lt_handle(t)+l+L_1-(L_1*L_2)/M;
        L_0=holderArray.Lt_handle(t);
        dL_dt=holderArray.dLt_handle(t);
        
        V_cap=y(2);
        I2=y(3);
        
        dI=zeros(3,1);
        dI(1,1)=-(I1*M*R_1 - I1*M*dL_dt + M*V_cap*a + I1*L_2*R_1*a^2 - ...
            I1*M*R_1*a^2 - I1*L_2*a^2*dL_dt + I1*M*a^2*dL_dt + I2*M*R_2*a)...
            /(M*l - M^2 - L_0*M + L_1*M - L_0*L_2*a^2 + L_1*L_2*a^2 + ...
            L_0*M*a^2 - L_1*M*a^2 + L_2*a^2*l - M*a^2*l);
        dI(2,1)=I2/C_load;
        dI(3,1)=-(a*(I1*M*R_1 - I1*M*dL_dt - L_0*V_cap*a + L_1*V_cap*a + ...
            V_cap*a*l - I2*L_0*R_2*a + I2*L_1*R_2*a + I2*R_2*a*l))...
            /(M*l - M^2 - L_0*M + L_1*M - L_0*L_2*a^2 + L_1*L_2*a^2 + ...
            L_0*M*a^2 - L_1*M*a^2 + L_2*a^2*l - M*a^2*l);
    end
end
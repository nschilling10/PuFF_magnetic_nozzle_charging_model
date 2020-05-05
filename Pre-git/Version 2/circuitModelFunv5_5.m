% Nathan Schilling
% 02/19/19
function [E_gain,E_circ] = circuitModelFunv5_5(circuitInputParams)
   
    % Physics constants
    mu0=4*pi*10^-7;
    
    % Nozzle Input parameters (SI units)
    x_0=1;
    A=pi*x_0^2;
    N=1e1;
    l_nozzle=1;
    v_0=(25e3/3);
    
    % Circuit parameters
    input.l_0=0;
    input.a=circuitInputParams.a;
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
    L0=circuitInputParams.L0;
    global tau
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

%     function out=L_nozz_hand(x)
%         s=0.5;
% %         if x < tau
%             out=L0*((-0.8/tau^2)*(1/(s^2-s))*x.^2+(0.8/tau)*(1/(s^2-s))*x+1);
% %         else
% %             out=0;
% %         end
%     end
% 
%     function out=dL_nozz_hand(x)
%         s=0.5;
% %         if x < tau
%             out=L0*(2*(-0.8/tau^2)*(1/(s^2-s))*x+(0.8/tau)*(1/(s^2-s)));
% %         else
% %             out=0;
% %         end
%     end

    % Set up input functions
    D2L_nozz_hand=@(x) 0*x;
    input.Lt_handle=@(x) L_nozz_hand(x);
    input.dLt_handle=@(x) dL_nozz_hand(x);
    
    % Inital Condition
    Iode=[0;I1_0];
    
    tSpan=linspace(0,tauPerc*tau,1e4);
    
    [t,y]=ode45(@(t,y) fcgfuns(t,y,input),tSpan,Iode);   

    I_1=y(:,2);
    V_Cap=y(:,1);
    
    I_2=input.a*I_1;
    
    voltage_in_circuit=L_nozz_hand(t).*[0;(diff(I_1)./diff(t))]+dL_nozz_hand(t).*I_1;
    
    if circuitInputParams.graphDisplay

        figure(1)
        semilogy(t*10^6,I_1,t*10^6,I_2)
        grid on
        xlabel('\textbf{Time since ignition, }\boldmath$\mu$\textbf{s}','interpreter','latex','fontsize',22)
        ylabel('textbf{Current, A}','interpreter','latex','fontsize',22)
        title('Current vs. Time since ignition with load connected')
        legend({'\boldmath $I_1$','\boldmath $I_2$'},'interpreter','latex','fontsize',18)
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
        
        E_gen=0.5*L_nozz_hand(t).*I_1.^2;
        E_l=0.5*input.l_0*I_1.^2;
        E_cap=0.5*input.C_load*V_Cap.^2;
        
        E_R1=ones(length(t),1);
        for i=1:length(E_R1)
            if i == 1
                E_R1(i)=0;
            else
                E_R1(i)=input.R1*trapz(t(1:i),I_1(1:i).^2);
            end
        end
        E_tot=E_gen+E_l+E_cap-E_R1;
    figure(5)
    plot(t*10^6,E_gen*10^-6,t*10^6,E_R1*10^-6,t*10^6,E_l*10^-6,t*10^6,E_cap*10^-6)
    xlabel('Time since ignition (\mus)')
    ylabel('Energy in component (MJ)')
    legend('Generator','R_1','l_0','Cap','Location','Northwest')
    figure(6)
    plot(t*10^6,E_tot*10^-6)
    xlabel('Time since ignition (\mus)')
    ylabel('Total energy in the circuit (MJ)')
    figure(7)
    plot(t*1e6,voltage_in_circuit)
    xlabel('Time since ignition (\mus)')
    ylabel('Voltage in circuit')
    end
    
    E_gain=0.5*input.C_load*(V_Cap(end)^2-V_Cap(1)^2);
    E_circ=0.5*L0*I_1(1)^2;
    
    %this is the ode circuit solver
    function dI = fcgfuns(t,y,holderArray)
        %varargin is a structure that is assumed to hold the following
        %fields
        V_cap=y(1);
        I1=y(2);
        
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
        
        if ~isfield(holderArray,'a')
            a=1;
        else
            a=holderArray.a;
        end
        
        if ~isfield(holderArray,'C_load')
            C_load=1;
        else
            C_load=holderArray.C_load;
        end
        
        L_0=holderArray.Lt_handle(t)-l;
        
        dI=zeros(2,1);
        dI(1,1)=I1/C_load;
        dI(2,1)=((R_1+a^2*R_2-holderArray.dLt_handle(t))*I1+a^2*V_cap)/L_0;
    end
end
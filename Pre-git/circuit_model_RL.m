function E_gain = circuit_model_RL(varargin)

    % Physics constants
    mu0=4*pi*10^-7;
    
    % Nozzle Input parameters (SI units)
    x_0=1;
    A=pi*x_0^2;
    N=1e1;
    l_nozzle=(3.5)^2/4;
    v_0=5e3;
    
    % Circuit parameters
    input.L1=1e-3;
    input.l_0=0;
    input.R1=0;
    
    % Seed current
    input.I1_0=5e3;
    
    % Display parameters
    tauPerc=1;
    
    % Nozzle inductive derived parameters
    global L0
    L0=mu0*N^2*A/l_nozzle;
    L0=0.2;
    global tau
    tau=x_0/v_0;
    tau=10e-6;
    
%     function out=L_nozz_hand(x)
%         if x < tau
%             out=L0*(1-x/tau);
%         else
%             out=0;
%         end
%     end
% 
%     function out=dL_nozz_hand(x)
%         if x < tau
%             out=L0*(-1/tau);
%         else
%             out=0;
%         end
%     end

    function out=L_nozz_hand(x)
        tau_1=2e-6;
        out=L0*(exp(-x.^2/(tau_1^2))+0.5*(tanh((x-tau_1*3)/tau_1*3)+1));
        
    end

    function out=dL_nozz_hand(x)
        tau_1=2e-6;
        out=-L0*((3*(tanh((3*(3*tau_1 - x))/tau_1).^2 - 1))/(2*tau_1) + ...
            (2*x.*exp(-x.^2/tau_1^2))/tau_1^2);
    end

    input.Lt_handle=@(x) L_nozz_hand(x);
    input.dLt_handle=@(x) dL_nozz_hand(x);
    
    % Inital Conditions
    Iode=[input.I1_0];
    
    tSpan=linspace(0,tauPerc*tau);
    
    [t,y]=ode45(@(t,y) fcgfuns(t,y,input),tSpan,Iode);   

    I_1=y(:,1);
    fontSize=36;
    figure(1)
    semilogy(t*10^6,I_1,'LineWidth',5)
    xlabel('\textbf{Time since ignition,} \boldmath$\mu$\textbf{s}')
    ylabel('Current (A)')
    title('Current vs. Time since ignition')
%    xlabel('\textbf{Time since ignition (}{\boldmath$\mu$}\textbf{sec)}','interpreter','latex','fontsize',fontSize)
%     ylabel('\textbf{Current (A)}','interpreter','latex','fontsize',fontSize)
%     set(gca,'fontsize',28)
%     set(gca,'LineWidth',2)
    grid on
    
    [E_gain,ind]=max(0.5*input.L1*I_1.^2);
    
    voltage_in_circuit=L_nozz_hand(t).*[0;(diff(I_1)./diff(t))]+dL_nozz_hand(t).*I_1;
    
    voltage_L1=input.L1.*[0;(diff(I_1)./diff(t))];
    
    figure(6)
    plot(t*1e6,0.5*input.L1*I_1.^2)
    xlabel('Time since ignition (\mus)')
    ylabel('Inductive energy')
    
    figure(8)
    plot(t*1e6,voltage_L1.*I_1)
    xlabel('Time since ignition (\mus)')
    ylabel('Power')
    Power_vec_L1=abs(voltage_L1.*I_1);
    E_heat_L1=zeros(length(Power_vec_L1),1);
    for i=2:length(Power_vec_L1)
        E_heat_L1(i)=trapz(t(1:i),Power_vec_L1(1:i));
    end
    figure(9)
    plot(t*1e6,E_heat_L1)
    xlabel('Time since ignition (\mus)')
    ylabel('Heating energy (J)')
    
    max_heat=max(E_heat_L1(ind))
    
    figure(7)
    plot(t*1e6,voltage_in_circuit)
    xlabel('Time since ignition (\mus)')
    ylabel('Voltage in circuit')
    
    %this is the ode circuit solver
    function dI = fcgfuns(t,y,holderArray)
        %varargin is a structure that is assumed to hold the following
        %fields
        I1=y;
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
            L_1=0;
        else
            L_1=holderArray.L1;
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
        
        if ~isfield(holderArray,'I1_0')
            I10=1e6;
        else
            I10=holderArray.I1_0;
        end

        L_0=holderArray.Lt_handle(t)-l-L_1;
        flux=I10*L0;
        dI = zeros(1,1);
        
        dI(1,1) = (1/L_0)*(R_1*I1-(holderArray.dLt_handle(t)*flux/holderArray.Lt_handle(t)));
    end
end
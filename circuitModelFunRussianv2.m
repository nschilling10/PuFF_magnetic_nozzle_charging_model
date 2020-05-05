% Nathan Schilling
% 06/03/19
function [E_gain,E_circ] = circuitModelFunRussianv2(circuitInputParams)
   
    % Physics constants
    mu0=4*pi*10^-7;
    
    % Nozzle Input parameters (SI units)
    x_0=1;
    A=pi*x_0^2;
    N=1e1;
    l_nozzle=1;
    v_0=(25e3/3);
    
    % Circuit parameters
    input.L1t=circuitInputParams.L1;
    input.L2t=circuitInputParams.Lratio*input.L1t;
    input.k=0.96;
    M_circ=input.k*sqrt(input.L1t*input.L2t);
    input.L_e=10e-9;
    input.R1=0;
    input.R2=input.R1;
    if ~isfield(circuitInputParams,'C_load')
        input.C_load=5e-6;
    else
        input.C_load=circuitInputParams.C_load;
    end
    
    % Display parameters
    tauPerc=1;
    
    % Seed current
    I1_0=0.5e6;
    
    % Nozzle inductive derived parameters
    global L0
    L0=1.4e-6;
    global tau
    tau=120e-6;
    
    function out=L_nozz_hand(x)
        if x < tau
            out=L0*(1-(Lg0/L0)*(x/tau));
        else
            out=0;
        end
    end

    function out=dL_nozz_hand(x)
        if x < tau
            out=Lg0*(-1/tau);
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
    
    tSpan=linspace(0,tauPerc*tau,1e4);
    
    [t,y]=ode45(@(t,y) fcgfuns(t,y,input),tSpan,Iode);   

    I_1=y(:,1);
    V_Cap=y(:,2);
    I_2=y(:,3);
    
    a=1./((tau*(L0/Lg0))-tSpan);
    b=1./(a*sqrt(input.L2t*input.C_load));
    %Z=(2*l/b*k^2)*sqrt((1+alpha)*(1+alpha-k^2*(1+a.*tSpan)/l));
    
    if circuitInputParams.graphDisplay

        figure(1)
        semilogy(t*10^6,I_1,t*10^6,I_2)
        grid on
        xlabel('Time since ignition (\mus)')
        ylabel('Current (A)')
        title('Current vs. Time since ignition with load connected')
        legend('Generator side current','Load side current','location','southeast')
        
%         fontSize=36;
%         xlabel('\textbf{Time since ignition (}{\boldmath$\mu$}\textbf{sec)}','interpreter','latex','fontsize',fontSize)
%         ylabel('\textbf{Current (A)}','interpreter','latex','fontsize',fontSize)
        
        figure(2)
        plot(t*10^6,I_1*1e-6)
        grid on
        xlabel('Time since ignition (\mus)')
        ylabel('Generator Current (MA)')
        title('Generator Current vs. Time since ignition')
%         xlabel('\textbf{Time since ignition (}{\boldmath$\mu$}\textbf{sec)}','interpreter','latex','fontsize',fontSize)
%         ylabel('\textbf{Current on generator side (MA)}','interpreter','latex','fontsize',fontSize)
%         set(gca,'fontsize',28)
%         set(gca,'LineWidth',2)
        
        figure(3)
        plot(t*10^6,y(:,2)*1e-3)
        grid on
        xlabel('Time since ignition (\mus)')
        ylabel('Capacitor voltage (kV)')
        title('Capacitor voltage vs. Time since ignition')
%         xlabel('\textbf{Time since ignition (}{\boldmath$\mu$}\textbf{sec)}','interpreter','latex','fontsize',fontSize)
%         ylabel('\textbf{Voltage across capacitor \textup{U}(kV)}','interpreter','latex','fontsize',fontSize)
%         set(gca,'fontsize',28)
%         set(gca,'LineWidth',2)

        figure(2)
        plot(t*10^6,I_2*1e-6)
        grid on
        xlabel('Time since ignition (\mus)')
        ylabel('Load Current (MA)')
        title('Load Current vs. Time since ignition')
%         xlabel('\textbf{Time since ignition (}{\boldmath$\mu$}\textbf{sec)}','interpreter','latex','fontsize',fontSize)
%         ylabel('\textbf{Current on generator side (MA)}','interpreter','latex','fontsize',fontSize)
%         set(gca,'fontsize',28)
%         set(gca,'LineWidth',2)
        
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
        
        if ~isfield(holderArray,'L1t')
            L_1t=1;
        else
            L_1t=holderArray.L1t;
        end
        
        if ~isfield(holderArray,'L2t')
            L_2t=1;
        else
            L_2t=holderArray.L2t;
        end
        
        if ~isfield(holderArray,'k')
            M=0.85*sqrt(L_1t*L_2t);
        else
            M=holderArray.k*sqrt(L_1t*L_2t);
        end
        
        if ~isfield(holderArray,'L_e')
            Le=0;
        else
            Le=holderArray.L_e;
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
        
        L1=holderArray.Lt_handle(t)+Le+L_1t;
        L2=L_2t;
        L12=holderArray.k*L_2t*sqrt(L_1t+Le);
        
        U=y(2);
        I2=y(3);
        
        dI=zeros(3,1);
        dI(1,1)=-(L12*U - I1*L2*R_1 + I2*L12*R_2 - I1*L2*holderArray.dLt_handle(t))/(L12^2 - L1*L2);
        dI(2,1)=I2/C_load;
        dI(3,1)=(L1*U + I2*L1*R_2 - I1*L12*R_1 - I1*L12*holderArray.dLt_handle(t))/(L12^2 - L1*L2);
%         dI(1,1)=(I1(1)*L2*R_1 - L12*U + I2*L12*R_2 + I1(1)*L2*holderArray.dLt_handle(t))/(L12^2 - L1*L2);
%         dI(2,1)=I2/C_load;
%         dI(3,1)=(I2*L1*R_2 - L1*U + I1(1)*L12*R_1 + I1(1)*L12*holderArray.dLt_handle(t))/(L12^2 - L1*L2);
    end
end
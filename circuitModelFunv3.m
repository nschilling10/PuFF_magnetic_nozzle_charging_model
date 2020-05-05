% Nathan Schilling
% 02/19/19
function [E_gain,E_circ] = circuitModelFunv3(circuitInputParams)
   
    % Physics constants
    mu0=4*pi*10^-7;
    
    % Nozzle Input parameters (SI units)
    x_0=1;
    A=pi*x_0^2;
    N=1e1;
    l_nozzle=1;
    v_0=5e3;
    
    % Circuit parameters
    input.L1=circuitInputParams.L1;
    input.L2=circuitInputParams.Lratio*input.L1;
    input.k=circuitInputParams.k;
    M_circ=input.k*sqrt(input.L1*input.L2);
    input.l_0=1e-9;
    input.R1=1e-3;
    input.R2=input.R1;
    input.C_load=(22e6/(0.5*(85e3)^2));
    
    % Display parameters
    tauPerc=0.90;
    
    % Seed current
    I1_0=5.7e5;
    %I1_0=1e6;
    
    % Nozzle inductive derived parameters
    global L0
    L0=mu0*N^2*A/l_nozzle;
    global tau
    tau=x_0/v_0;
    
    function out=L_nozz_hand(x)
        if x < tau
            out=L0*(1-x/tau);
        else
            out=0;
        end
    end

    function out=dL_nozz_hand(x)
        if x < tau
            out=L0*(-1/tau);
        else
            out=0;
        end
    end

    % Set up input functions
    D2L_nozz_hand=@(x) 0*x;
    input.Lt_handle=@(x) L_nozz_hand(x);
    input.dLt_handle=@(x) dL_nozz_hand(x);
    
    % Inital Conditions??
    Iode=[I1_0;0;0];
    
    tSpan=linspace(0,tauPerc*tau,1e4);
    %tSpan=[linspace(0,1,1e5)];
    %tSpan=[tSpan linspace(1+tSpan(2)-tSpan(1),1e5,1e5)];
    
    [t,y]=ode45(@(t,y) fcgfuns(t,y,input),tSpan,Iode);   

    I_1=y(:,1);
    V_Cap=y(:,2);
    I_2=y(:,3);
    
    if circuitInputParams.graphDisplay
        fontSize=14;
        figure(1)
        semilogy(t,I_1,t*10^6,I_2,'LineWidth',5)
        grid on
        xlabel('\textbf{Time since ignition (}\textbf{sec)}','interpreter','latex','fontsize',fontSize)
        ylabel('\textbf{Current (A)}','interpreter','latex','fontsize',fontSize)
        set(gca,'fontsize',fontSize+4)
        set(gca,'LineWidth',2)
        
        figure(2)
        plot(t,I_1*1e-6,'LineWidth',5)
        grid on
        xlabel('\textbf{Time since ignition (}\textbf{sec)}','interpreter','latex','fontsize',fontSize)
        ylabel('\textbf{Current on generator side (MA)}','interpreter','latex','fontsize',fontSize)
        set(gca,'fontsize',fontSize+4)
        set(gca,'LineWidth',2)
        
        figure(3)
        plot(t,y(:,2)*1e-3,'LineWidth',5)
        grid on
        xlabel('\textbf{Time since ignition (}\textbf{sec)}','interpreter','latex','fontsize',fontSize)
        ylabel('\textbf{Voltage across capacitor \textup{U}(kV)}','interpreter','latex','fontsize',fontSize)
        set(gca,'fontsize',fontSize+4)
        set(gca,'LineWidth',2)
    end
    
    E_gain=0.5*input.C_load*(V_Cap(end)^2-V_Cap(1)^2);
    E_circ=0.5*L0*I_1(1)^2;
    
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
        L_0=l+L_1-holderArray.Lt_handle(t);
        
        V_cap=I1(2);
        I2=I1(3);
        
        dI=zeros(3,1);
        dI(1,1)=(M*V_cap + I1(1)*L_2*R_1 + I2*M*R_2 - I1(1)*L_2*holderArray.dLt_handle(t))/(M^2 - L_0*L_2);
        dI(2,1)=I2/C_load;
        dI(3,1)=(L_0*V_cap + I2*L_0*R_2 + I1(1)*M*R_1 - I1(1)*M*holderArray.dLt_handle(t))/(M^2 - L_0*L_2);
    end
end
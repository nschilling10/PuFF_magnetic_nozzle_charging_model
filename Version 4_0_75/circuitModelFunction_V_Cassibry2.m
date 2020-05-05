% Nathan Schilling
% 03/31/20
% outputs t,I_1,I_2,V_Cap,Rp_vec, Vp_vec
% Requires L1, L2 (transformer inductances) R1, R2, magnetic coupling
% constant k, C, R_fcc (radius of the nozzle), inital current, V_sqig, and P_mag. 
% V_sqig=dB/dt * A (no dL/dt or L(t))
function [t,I_1_vec,I_2_vec,V_Cap,Rp_vec,Vp_vec] = circuitModelFunction_V_Cassibry2(circuitInputParams)
    
    % Circuit parameters
    input.L_T1=circuitInputParams.L1;
    input.L_T2=circuitInputParams.L2;
    input.k=circuitInputParams.k;
    M_circ=input.k*sqrt(input.L_T1*input.L_T2);
    input.L_1=circuitInputParams.l_1;
    input.L_2=circuitInputParams.l_2;
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
    input.C_load=circuitInputParams.C;
    input.R_Fcc=circuitInputParams.R_Fcc;
    input.P_mag_hand=@(r_p) circuitInputParams.P_mag(r_p);
    input.m_p=circuitInputParams.m_p;
    
    % Display parameters
    tauPerc=circuitInputParams.tauPerc;
    
    % Inital values
    I1_0=circuitInputParams.I0;
    Rp_0=circuitInputParams.Rp0;
    Vp_0=circuitInputParams.vp0;
    
    % Set up input functions
    tau=circuitInputParams.tau;
    input.V_sqig_hand=circuitInputParams.V_sqig;
    
    % Inital Condition [I1_0, Vcap_0, I2_0, Rp_0, vp_0]
    Iode=[I1_0;0;0;Rp_0;Vp_0];
    
    tSpan=linspace(0,tauPerc*tau,1e3);
    
    options=odeset('Events',@zeroRpStopEvent);
    [t,y]=ode45(@(t,y) fcgfuns(t,y,input),tSpan,Iode);   
    
%     %options=odeset('RelTol',1e-20,'AbsTol',1e-20);
%     options=odeset('Events',@capBackStopEvent);
%     
%     soln=ode45(@(t,y) fcgfuns(t,y,input),tSpan,Iode,options);   
%     y=soln.y';
%     t=soln.x;
%     
%     % Run again with R_2=1M Ohm
%     input.R2=1e6;
%     options=odeset('Events',@zeroRpStopEvent);
%     soln=ode45(@(t,y) fcgfuns(t,y,input),linspace(t(end),tauPerc*tau,1e3),y(end,:)',options);
%     y=[y;soln.y'];
%     t=[t soln.x];
    
    I_1_vec=y(:,1);
    V_Cap=y(:,2);
    I_2_vec=y(:,3);
    Rp_vec=y(:,4);
    Vp_vec=y(:,5);
    
    %this is the ode circuit solver
    function dI = fcgfuns(t,y,holderArray)
        %varargin is a structure that is assumed to hold the following
        %fields
        I_1=y(1);
        V_cap=y(2);
        I_2=y(3);
        r_p=y(4);
        v_p=y(5);
        
        if ~isfield(holderArray,'V_sqig_hand')
            msgID = 'fcgfuns:BadInput';
            msg = 'No V_sqig_hand function specified.';
            baseException = MException(msgID,msg);
            throw(baseException);
        else
            V_sqig=holderArray.V_sqig_hand(r_p,v_p);
        end
        
        if ~isfield(holderArray,'P_mag_hand')
            msgID = 'fcgfuns:BadInput';
            msg = 'No P_mag function specified.';
            baseException = MException(msgID,msg);
            throw(baseException);
        else
            %P_mag=holderArray.P_mag_hand(I1);
            P_mag=holderArray.P_mag_hand(r_p);
        end
        
        if ~isfield(holderArray,'L_T1')
            L_T1=1;
        else
            L_T1=holderArray.L_T1;
        end
        
        if ~isfield(holderArray,'L_T2')
            L_T2=1;
        else
            L_T2=holderArray.L_T2;
        end
        
        if ~isfield(holderArray,'k')
            M=0.9*sqrt(L_T1*L_T2);
        else
            M=holderArray.k*sqrt(L_T1*L_T2);
        end
        
        if ~isfield(holderArray,'L_1')
            L_1=0;
        else
            L_1=holderArray.L_1;
        end
        
        if ~isfield(holderArray,'L_2')
            L_2=0;
        else
            L_2=holderArray.L_2;
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
        
        if ~isfield(holderArray,'m_p')
            m_p=0;
        else
            m_p=holderArray.m_p;
        end        
                
        if ~isfield(holderArray,'R_Fcc')
            R_Fcc=0;
        else
            R_Fcc=holderArray.R_Fcc;
        end
        
        dI=zeros(5,1);
        dI(1,1)=-(M*V_cap - L_T2*V_sqig - L_2*V_sqig + I_1*L_2*R_1 + I_1*L_T2*R_1 + I_2*M*R_2)...
            /(- M^2 + L_1*L_2 + L_1*L_T2 + L_2*L_T1 + L_T1*L_T2);
        dI(2,1)=I_2/C_load;
        dI(3,1)=-(L_1*V_cap + L_T1*V_cap - M*V_sqig + I_2*L_1*R_2 + I_2*L_T1*R_2 + I_1*M*R_1)...
            /(- M^2 + L_1*L_2 + L_1*L_T2 + L_2*L_T1 + L_T1*L_T2);
        dI(4,1)=v_p;
        dI(5,1)=(-P_mag/m_p)*(2*pi*r_p*R_Fcc);
    end

    % Stops integration if Rp is 0
    function [value,isterminal,direction] = zeroRpStopEvent(t,y)
        value=y(4);
        isterminal=1;
        direction=[];
    end
    
    % Stop the integration if the voltage on the capacitor starts to
    % decrease
    function [value,isterminal,direction] = capBackStopEvent(t,y)
        value=y(3);
        isterminal=1;
        direction=-1;
    end
end
% Nathan Schilling
% 11/24/19
function [t,I_1,I_2,V_Cap] = circuitModelFunction_V_Cassibry(circuitInputParams)
    
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
    
    % Display parameters
    tauPerc=1;
    
    % Seed current
    I1_0=circuitInputParams.I0;
    
    % Set up input functions
    tau=circuitInputParams.tau;
    input.v_handle=@(t) circuitInputParams.v_exp_hand(t);
    input.Lt_handle=@(d) circuitInputParams.L_nozz_hand(d);
    input.dLt_handle=@(t) circuitInputParams.dL_nozz_hand(t);
    
    % Inital Condition [I1_0, Vcap_0, I2_0, Rp_0]
    Iode=[I1_0;0;0;circuitInputParams.Rp0];
    
    tSpan=linspace(0,tauPerc*tau,1e3);
    
    options=odeset('RelTol',1e-20,'AbsTol',1e-20);
    
    [t,y]=ode45(@(t,y) fcgfuns(t,y,input),tSpan,Iode);   

    I_1=y(:,1);
    V_Cap=y(:,2);
    I_2=y(:,3);
    Rp=y(:,4);
    d_vec=0.5*holderArray.R_Fcc-Rp;
    L_t_vec=holderArray.Lt_handle(d_vec);
    
    figure
    plot(tSpan,L_t_vec)
    keyboard
    
    %this is the ode circuit solver
    function dI = fcgfuns(t,y,holderArray)
        %varargin is a structure that is assumed to hold the following
        %fields
        I1=y(1);
        V_cap=y(2);
        I2=y(3);
        d=0.5*holderArray.R_Fcc-y(4);
        
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
        
        if ~isfield(holderArray,'v_handle')
            msgID = 'fcgfuns:BadInput';
            msg = 'No dLt_handle input field found.';
            baseException = MException(msgID,msg);
            throw(baseException);
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
        
        V_sqig=I1*holderArray.dLt_handle(t);
        
        L_0=L_T1+holderArray.Lt_handle(d)+L_1;
        
        di2_dt=(V_cap-I2*R_2-(M*(V_sqig-I1*R_1)/L_0))/(L_2+L_T2-(M^2/L_0));
        
        dI=zeros(3,1);
        dI(1,1)=(V_sqig-I1*R_1-M*di2_dt)/(L_0);
        dI(2,1)=I2/C_load;
        dI(3,1)=di2_dt;
        dI(4,1)=holderArray.v_handle(t);
    end
end
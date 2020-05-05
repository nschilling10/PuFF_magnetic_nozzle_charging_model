% Nathan Schilling
% 02/19/19
function [t,I_1,I_2,V_Cap,d_vec] = circuitModelFunValidatedv2_1(circuitInputParams)
   
    % Circuit parameters
    input.L1=circuitInputParams.L1;
    input.L2=circuitInputParams.L2;
    input.k=circuitInputParams.k;
    M_circ=input.k*sqrt(input.L1*input.L2);
    input.l_1=circuitInputParams.l_1;
    input.l_2=circuitInputParams.l_2;
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
    
    tau=circuitInputParams.tau;
    
    % Set up input functions
    input.Lt_handle=circuitInputParams.L_nozz_hand;
    input.dLt_handle=circuitInputParams.dL_nozz_hand;
    input.v_handle=circuitInputParams.v_exp_hand;
    
    % Inital Condition [I1_0, V_cap_0, I2_0, Rp_0]
    Iode=[I1_0;0;0;circuitInputParams.Rp0];
    
    tSpan=linspace(0,tauPerc*tau,1e3);
    
    [t,y]=ode45(@(t,y) fcgfuns(t,y,input),tSpan,Iode);
    
    I_1=y(:,1);
    V_Cap=y(:,2);
    I_2=y(:,3);
    Rp=y(:,4);
    d_vec=0.5*input.R_Fcc-Rp;
    
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
        
        if ~isfield(holderArray,'l_1')
            l_1=0;
        else
            l_1=holderArray.l_1;
        end
        
        if ~isfield(holderArray,'l_2')
            l_2=0;
        else
            l_2=holderArray.l_2;
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
        
        L_01=L_1+l_1+holderArray.Lt_handle(d);
        L_02=L_2+l_2;
        
        dI=zeros(4,1);
        
        dI(1,1)=(M*V_cap + I1*L_02*R_1 + I2*M*R_2 + I1*L_02*holderArray.dLt_handle(t))/(M^2 - L_01*L_02);
        dI(2,1)=I2/C_load;
        dI(3,1)=(L_01*V_cap + I2*L_01*R_2 + I1*M*R_1 + I1*M*holderArray.dLt_handle(t))/(M^2 - L_01*L_02);
        dI(4,1)=holderArray.v_handle(t);
    end
end
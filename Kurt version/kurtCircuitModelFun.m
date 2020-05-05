function [E_gen,E_cap]=kurtCircuitModelFun(circuitInputParams)
    
    % Nozzle inductive derived parameters
    global L0
    %L0=mu0*N^2*A/l_nozzle;
    L0=circuitInputParams.L0;
    global tau
    %tau=x_0/v_0;
    tau=120e-6;

%     function out=L_nozz_hand(x)
%         tau_1=2e-6;
%         %         if x < tau
%         out=L0*(exp(-x.^2/(tau_1^2))+0.5*(tanh((x-tau_1*3)/tau_1*3)+1));
%         %         else
%         %             out=0;
%         %         end
%     end
% 
%     function out=dL_nozz_hand(x)
%         tau_1=2e-6;
%         %         if x < tau
%         out=-L0*((3*(tanh((3*(3*tau_1 - x))/tau_1).^2 - 1))/(2*tau_1) + ...
%             (2*x.*exp(-x.^2/tau_1^2))/tau_1^2);
%         %         else
%         %             out=0;
%         %         end
%     end


    struct4Input.L_gen=@(x) L_nozz_hand(x);
    struct4Input.dL_gen=@(x) dL_nozz_hand(x);

    struct4Input.Lp=circuitInputParams.Lp;
    struct4Input.Lc=circuitInputParams.Lc;
    struct4Input.L0=0;
    struct4Input.Rp=circuitInputParams.Rp;
    struct4Input.Rc=circuitInputParams.Rc;
    struct4Input.C_load=circuitInputParams.C_load;
    
    struct4Input.k=0.85;
    
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

    tSpan=linspace(0,tau,1e4);
    
    I0=[0;circuitInputParams.I0;0];
    
    [t,y]=ode45(@(t,y) f(t,y,struct4Input),tSpan,I0);
    
    I1_vec=y(:,1);
    I2_vec=y(:,2);
    V_cap_vec=y(:,3);
    
    figure(1)
    plot(t*1e6,I1_vec*1e-6)
    xlabel('t in microseconds')
    ylabel('I_1 (MA)')
    
    figure(2)
    plot(t*1e6,I2_vec*1e-6)
    xlabel('t in microseconds')
    ylabel('I_2 (MA)')
    
    figure(3)
    plot(t*1e6,V_cap_vec*1e-3)
    xlabel('t in microseconds')
    ylabel('V_cap (kV)')
    
    E_gen=0.5*L0*I2_vec(1)^2;
    E_cap=0.5*circuitInputParams.C_load*V_cap_vec(end)^2;
    
    voltage_in_circuit=L_nozz_hand(t).*[(diff(I2_vec)./diff(t));0]+dL_nozz_hand(t).*I2_vec;
    
    v_Lc=struct4Input.Lc.*[(diff(I1_vec)./diff(t));0];
    v_Lp=struct4Input.Lp.*[(diff(I2_vec)./diff(t));0];
    v_L0=struct4Input.L0.*[(diff(I1_vec)./diff(t));0];
    v_Rc=struct4Input.Rc*I1_vec;
    v_Rp=struct4Input.Rp*I2_vec;
    M_circ=struct4Input.k*sqrt(struct4Input.Lp*struct4Input.Lc);
        
    V_loop1= v_Lc+v_Rc - M_circ*([(diff(I2_vec)./diff(t));0]) + v_L0 + V_cap_vec;
    V_loop2= v_Lp+v_Rp - M_circ*([(diff(I1_vec)./diff(t));0]) + voltage_in_circuit;

    figure(7)
    plot(t*1e6,V_loop1,t*1e6,V_loop2)
    xlabel('Time since ignition (\mus)')
    ylabel('Voltage in a loop')
    legend('V loop 1','V loop 2')
   
    function dI = f(t,y,inputStruct)
       L_p=inputStruct.Lp;
       L_c=inputStruct.Lc;
       L_0=inputStruct.L0;
       R_p=inputStruct.Rp;
       R_c=inputStruct.Rc;
       C=inputStruct.C_load;
       M=inputStruct.k*sqrt(L_p*L_c);
       
       I1=y(1);
       I2=y(2);
       V=y(3);
       
       L_p=L_p+inputStruct.L_gen(t);
       R_p=R_p+inputStruct.dL_gen(t);
       
       dI(1,1)=(V*L_p-I2*M*R_p-I1*R_c*L_p)/((L_p+L_0)*L_c - M^2);
       dI(2,1)=(M*dI(1,1)-I2*R_p)/(L_p);
       dI(3,1)= - I2/C;
    end
end
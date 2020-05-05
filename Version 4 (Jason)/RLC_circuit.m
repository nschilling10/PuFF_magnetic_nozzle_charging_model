% Nathan Schilling
% 11/24/19
function RLC_circuit()
    
    % Circuit parameters
    input.R_load=5e-3;
    input.R_loss=0;
    input.R=input.R_load+input.R_loss;
    
    input.L_load=12e-9;
    input.L_loss=0;
    input.L=input.L_load+input.L_loss;
    
    input.C=400e-6;
    
    tau=(2*pi*sqrt(input.L*input.C));
    tauPerc=0.5;
    
    V0=100e3;
    
    % Inital Condition [I1_0, Vcap_0]
    Iode=[0;V0];
    
    tSpan=linspace(0,tauPerc*tau,1e3);
    
    options=odeset('RelTol',1e-10,'AbsTol',1e-10,'Events',@eventsFcn);
    
    [t,y]=ode45(@(t,y) fcgfuns(t,y,input),tSpan,Iode);   

    I_1=y(:,1);
    V_Cap=y(:,2);
    
    figure(1)
    plot(t*10^6,I_1*1e-6)
    grid on
    xlabel('\textbf{Time, }\boldmath$\mu$\textbf{s}','interpreter','latex','fontsize',22)
    ylabel('\textbf{Current, MA}','interpreter','latex','fontsize',22)
    title('\textbf{Current vs. Time since load connected}','interpreter','latex','fontsize',22)
    h.Children(2).LineWidth=2;
    h.Children(2).FontSize=18;
    set(findall(gca, 'Type', 'Line'),'LineWidth',2);
    
    figure(2)
    plot(t*10^6,V_Cap*1e-3)
    grid on
    title('\textbf{U(kV) vs. Time since ignition}','interpreter','latex')
    xlabel('\textbf{Time since ignition (}{\boldmath$\mu$}\textbf{sec)}','interpreter','latex')
    ylabel('\textbf{Voltage across capacitor \textup{U}(kV)}','interpreter','latex')
    set(gca,'fontsize',16)
    set(findall(gca, 'Type', 'Line'),'LineWidth',2);
    ylim([V_Cap(end)*1e-3 V_Cap(1)*1e-3])
    
    %this is the ode circuit solver
    function dI = fcgfuns(t,y,holderArray)
        I1=y(1);
        V_cap=y(2);
        
        dI=zeros(2,1);
        dI(1,1)=(V_cap-I1*(holderArray.R))/holderArray.L;
        dI(2,1)=-I1/holderArray.C;
    end

    %this is an event function
    function [position,isterminal,direction] = eventsFcn(t,y)
        position = y(2); % The value that we want to be zero
        isterminal = 1;  % Halt integration
        direction = -1;   % The zero can be approached from either direction
    end
end
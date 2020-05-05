function rlc_recharge_nozzle_circuit_v2()

%units
eV = 11605;
keV = 1000*eV;

%constants
mu0 = 4*pi*1e-7;


%inputs
m_propellant = .01;

if 0 %case used in lecture notes, 40 kV charged to caps
    %primary circuit (1) inputs
    L1 = 100e-9;
    R1 = 5e-3;
    
    %secondary circuit (2) inputs including reactor capacitance C
    L2 = 100e-9;
    R2 = 5e-3;
    C = 100e-6;
    
    %plasma inputs
    MW = 4;
    Rgas = 8314.5/MW;
    T0 = .5*keV;
    gamma = 1.3;
    
    %transformer inputs
    mu_r = 1;
    
    rT = .1;
    AT = rT^2 * pi;
    NT1 = 10;
    NT2 = NT1/10;
    lT1 = rT*1;
    lT2 = rT*3;
    
    %flux compression generator inputs
    Nfcc = 10;
    Rfcc = .95; %1.1;  %coil radius
    taufcc = Rfcc/(2*sqrt(2*gamma*Rgas*2*T0));
    
    %calculations from inputs
    LT1 = mu0 * mu_r * AT*NT1^2/lT1;
    LT2 = mu0 * mu_r * AT*NT2^2/lT2;
    MT = .9*sqrt(LT1*LT2);   
    
    %-----------initial conditions
    I10 = .51e6; %Amps, initial current in primary
    I20 = 0;
%     V10 = 0;
    V20 = 0;  %voltage on reactor caps
elseif 1  %PUFF
    m_propellant = .125;
    %primary circuit (1) inputs
    L1 = 100e-9;
    R1 = 0*1e-5 ;
    
    %secondary circuit (2) inputs including reactor capacitance C
    L2 = 10e-9;
    R2 = 5e-3 * 1;
    C = .01;  %because we need 50 MJ or more
    
    %plasma inputs
    MW = 4;
    Rgas = 8314.5/MW;
    T0 = 1*keV*.98;
    gamma = 1.3;
    
    %transformer inputs
    mu_r = 1;
    
    rT = .1;
    AT = rT^2 * pi;
    NT1 = 30;
    NT2 = 1.7; %NT1/10;
    lT1 = rT*1;
    lT2 = rT*4;
    
    %flux compression generator inputs
    Nfcc = 9;
    Rfcc = 18;
    taufcc = Rfcc/(2*sqrt(2*gamma*Rgas*2*T0));
    
    %calculations from inputs
    LT1 = mu0 * mu_r * AT*NT1^2/lT1;
    LT2 = mu0 * mu_r * AT*NT2^2/lT2;
    k=0.9;
    MT = k*sqrt(LT1*LT2);
        
    
    %-----------initial conditions
    I10 = .75e6; %Amps, initial current in primary
    I20 = 0;
%     V10 = 0;
    V20 = 0;  %voltage on reactor caps
else  %reactor grade caps
    
    m_propellant = .01;
    %primary circuit (1) inputs
    L1 = 100e-9;
    R1 = 5e-3 ;
    
    %secondary circuit (2) inputs including reactor capacitance C
    L2 = 100e-9;
    R2 = 5e-3 * 1;
    C = 400e-6;
    
    %plasma inputs
    MW = 4;
    Rgas = 8314.5/MW;
    T0 = 1*keV;
    gamma = 1.3;
    
    %transformer inputs
    mu_r = 1;
    
    rT = .1;
    AT = rT^2 * pi;
    NT1 = 20;
    NT2 = 2;
    lT1 = rT*1;
    lT2 = rT*3;
    
    %flux compression generator inputs
    Nfcc = 10;
    Rfcc = 3;
    taufcc = Rfcc/(2*sqrt(2*gamma*Rgas*2*T0));
    
    %calculations from inputs
    LT1 = mu0 * mu_r * AT*NT1^2/lT1;
    LT2 = mu0 * mu_r * AT*NT2^2/lT2;
    MT = .9*sqrt(LT1*LT2);
        
    
    %-----------initial conditions
    I10 = .77e6; %Amps, initial current in primary
    I20 = 0;
%     V10 = 0;
    V20 = 0;  %voltage on reactor caps
end

%make sure winding is backward.
% MT = -abs(MT);

%calculation of secondary initial conditions
B0 = mu0 * Nfcc * I10 /Rfcc;  %initial magnetic field, this is an approximation using the radius of the coils as a length scale
Phi0 = B0 * pi * Rfcc.^2;  %assuming this is conserved
A0 = pi*Rfcc.^2; %cross sectional area of nozzle


%% functions

% flux compression voltage model
% vexpfun = @(t) sqrt(gamma*Rgas*2*T0).* (cos(pi/2*t/taufcc) .* ((t/taufcc)<2)  -  (1.*(t/taufcc)>=2)) ;
Vfccfun = @(I,v) mu0 .* I * Nfcc^2/2 .* v;

%gap between field coil and plasma surface
dfun = @(r) Rfcc/2 - r;

%% initial conditions
Rp0 = 0;
vexp0 = sqrt(gamma*Rgas*2*T0);

%voltage on coil from expanding plasma
V10 = Vfccfun(I10,vexp0);


%time constant, based on flux compression time
tf = taufcc*20; %   Rfcc/vexp0 * 3;
dt = 100e-9;
%initial conditions
%                     initial plasma radius, initial plasma expansion rate
IV = [I10 I20 V10 V20 Rp0                    vexp0];

%options to ode
options = odeset('Events',@propellant_escapes_EventsFcn,'MaxStep',10e-9);

%solution to ode
[t,y,te,ye,ie]=ode45(@odefuns,[0:dt:tf],IV,options);



%call ode again if needed
if ie(end)==2 %event triggered by plasma hitting wall or being redirected by other nozzle field, not flux coil field, need to redirect flow
    tout = t; yout = y; teout = te; yeout = ye; ieout = ie;

    %initial conditions, which need some modifications
    IV2 = y(end,:);  
    %redirect the propellant flow
    IV2(6) = -IV2(6);
    %voltage on coil from expanding plasma
    IV2(3) = Vfccfun(IV2(1),vexp0);

    disp('Plasma bounced at nozzle wall');
    [t,y,te,ye,ie]=ode45(@odefuns,[t(end) tf],IV2,options);

    % Accumulate output.  
    nt = length(t);
    t = [tout; t(2:nt)];
    y = [yout; y(2:nt,:)];
    te = [teout; te];          % Events at tstart are never reported.
    ye = [yeout; ye];
    ie = [ieout; ie];

end


%call ode again if needed
if any(ie(2:end)==3) %event triggered by capacitor reaching its peak voltagae
    tout = t; yout = y; teout = te; yeout = ye; ieout = ie;

    %initial conditions, which need some modifications
    IV2 = y(end,:);  
    
    %turn off transformer and secondary circuit
    LT1 = 0;
    LT2 = 0;
    MT = 0;
    R2 = 1e0;
    
    disp('Capacitor reached peak charging voltage');
    [t,y,te,ye,ie]=ode45(@odefuns,[t(end) tf],IV2,options);

    % Accumulate output.  
    nt = length(t);
    t = [tout; t(2:nt)];
    y = [yout; y(2:nt,:)];
    te = [teout; te];          % Events at tstart are never reported.
    ye = [yeout; ye];
    ie = [ieout; ie];
    
end




I1 = y(:,1);
I2 = y(:,2);
R_plasma = y(:,5);
V_plasma = y(:,6);

% Vfcc = y(:,3);
Vfcc = Vfccfun(I1,V_plasma);
Vc = y(:,4);


% ---------- Force on plasma vs time
%area between coil and plasma
Area = A0 - 2*pi*R_plasma.^2;  %I'm assuming the plasma is a cylinder equal to Rp, it's current radius

%updated magnetic field between plasma
Bfield = Phi0./Area;

%calculate the force and acceleration acting on the plasma
Force = -Bfield.^2 /2/mu0 .* Area;

% ---------- Energies in the circuit
MJ = 1e-6;
% E_fcc =  .5 *mu0 * Nfcc^2  * dfun(R_plasma)/2 .* I1.^2 *MJ;
E_fcc =  .5 *mu0 * Nfcc^2  * Rfcc/2 .* I1.^2 *MJ;
E_plasma = .5*m_propellant*V_plasma.^2 *MJ;
E_cap = .5*C*Vc.^2 *MJ;
E_tot = E_fcc + E_plasma + E_cap;

Gain=E_cap(end)/(E_plasma(end)-E_plasma(1))

%plot the results
figure(1),
plot(t*1e6,I1/1e6), hold on
plot(t*1e6,I2/1e6), hold on
plot(t*1e6,Vfcc/1e6), hold off
xlabel('t (\mus)');
ylabel('I (MA), V(MV)')
lh = legend('$$I_1$$','$$I_2$$','$$\tilde{V}$$','location','southeast'); legend boxoff;
lh.Interpreter = 'latex';
grid on

figure(2),
plot(t*1e6,Vc/1e3)
xlabel('t (\mus)');
ylabel('V_c (kV)')
grid on
% lh = legend('Recharge Circuit (kV)','Flux Compression (MV)'); legend boxoff;

figure(3),
plot(t*1e6,R_plasma/(Rfcc/2))
xlabel('t (\mus)');
ylabel('R_{plasma}/R_{coil}')
grid on

figure(4),
plot(t*1e6,V_plasma/1e3)
xlabel('t (\mus)');
ylabel('V_{plasma} (km/s)')
grid on

figure(5),
plot(t*1e6,Bfield)
xlabel('t (\mus)');
ylabel('B (T)')
grid on

%energy in plasma and circuit
figure(6),
semilogy(t*1e6,E_fcc, ...
     t*1e6,E_plasma, ...
     t*1e6,E_cap, ...
     t*1e6,E_tot)
xlabel('t (\mus)');
ylabel('E (MJ)')
grid on
lh6 = legend('Coil','Plasma','Cap','Total'); legend boxoff
ylim([1e-3 1e6])

% keyboard

% figure(3)
% plot(t/taufcc,sqrt(gamma*Rgas*2*T0) * cos(2*pi*t/4/taufcc)/1000);
% 
% keyboard

% [ax,h1,h2] = plotyy(t*1e6,Vc/1e3,t*1e6,Vfcc/1e6); %, hold on 
% % plot(t*1e6,Vc/1e3), hold off
% xlabel('t (\mus)');
% ylabel(ax(1),'Capacitor Voltage (kV)')
% ylabel(ax(2),'Flux Compression Voltage (MV)')


%keyboard4


    function dy = odefuns(t,y)
        %
        %y1 = I1
        %y2 = I2
        %y3 = V1
        %y4 = V2
        %y5 = plasma radius
        %y6 = plasma velocity

        dy = zeros(length(IV),1);
 
        %-------------update variables
        I1 = y(1);
        I2 = y(2);
        V2 = y(4);
        %         V1 = -(mu0 * I1 * Nfcc * sqrt(2*gamma*Rgas*T0))/2 *cos(pi*t/2/taufcc);
        %expansion velocity, plasma/coil gap distance d
%         vexp = vexpfun(t);
        vexp = y(6); %expansion velocity

        
        %plasma geometry
        Rp = max(y(5),0);  %plasma radius
        d = dfun(Rp);
                        
        %area between coil and plasma
        A = A0 - 2*pi*Rp*Rp;  %I'm assuming the plasma is a cylinder equal to Rp, it's current radius

        %update reference field and flux from coil based on updated coil current, without plasma
        Bref = mu0 * Nfcc * I1 /Rfcc;  %initial magnetic field, this is an approximation using the radius of the coils as a length scale
        Phi = Bref * pi * Rfcc.^2;  %assuming this is conserved
        A0 = pi*Rfcc.^2; %max cross sectional area of nozzle
        
        %updated magnetic field between plasma
        B = Phi/A;
                
        %calculate the force and acceleration acting on the plasma
        F = -B^2 /2/mu0 * A;
        
        %turn on MHD, which acts like a drag term on the plasma and a voltage term on the primary coil
        if vexp<0 %plasma is rebounding and going out the exhaust.  Turn on the MHD generator.  Whoopee doo!
            VMHD = abs(vexp) * 1 * 2*pi*Rfcc * (I1<I10);
            dudt_mhd = -abs(I1*VMHD/m_propellant./(abs(vexp) + .1) *(abs(vexp)>100))* (I1<I10);
        else
            VMHD = 0;
            dudt_mhd = 0;
        end
        
        %
        %         V1 = y(3); %Vfccfun(I1,t,vexp);
        V1 = Vfccfun(I1,vexp) + VMHD;
        dudt = F/m_propellant*(Rp>0) + dudt_mhd;
            
        %flux compression coil fcc inductance
        Lfcc = mu0 * Nfcc^2 *d/2;
        
        %circuit equations
        
        %voltage on secondary
        dy(4) = -I2/C;
        
        %current on secondary
        dy(2) = (V2 - I2 * R2 - MT *((V1 - I1*R1)/(L1+Lfcc+LT1)))/(L2 + LT2 - (MT^2/(L1+Lfcc+LT1)) );
                
        
        %voltage on primary, which I can't make work right as an ode because you would have to include time derivative of I1 and have a 2nd term below for dy(3)
        dy(3) = I1*mu0*Nfcc^2 * (-dudt)*0;
        
        %current on primary
        dy(1) = (V1-I1*R1-MT*dy(2))/(L1+Lfcc+LT1);
        
        %plasma radius and expansion rate
        dy(5) = vexp;  
        dy(6) = dudt;
        
%         if t>(taufcc*.95)
%             keyboard
%         end
%                 
    end

    function [position,isterminal,direction] = propellant_escapes_EventsFcn(t,y)
        position(1) = y(5); % plasma radius, The value that we want to be zero
        isterminal(1) = 1;  % Halt integration if plasma leaves nozzle
        direction(1) = -1;   % The zero can be approached from negative direction (after flow has been redirected)

        position(2) =  dfun(y(5)) - .05*dfun(0); % The value that we want to be zero, when gap is less than 5% of original value
        isterminal(2) = 1;  % Halt integration if plasma hits the wall protecting the coils and is abruptly redirected
        direction(2) = 0;   % The zero can be approached from either direction
        
        %capacitor charged up as much as it will get
        position(3) =  y(2); % The value that we want to be zero, when gap is less than 5% of original value
        isterminal(3) = 2;  % Halt integration if plasma hits the wall protecting the coils and is abruptly redirected
        direction(3) = 0;   % The zero can be approached from either direction

        %flux coil current falls below seed current value
        position(4) =  y(1)-I10; % The value that we want to be zero, when gap is less than 5% of original value
        isterminal(4) = 1;  % Halt integration if plasma hits the wall protecting the coils and is abruptly redirected
        direction(4) = -1;   % The zero can be approached only when the coil current is decreasing from its peak value

        
        
        %         if t>1e-6
%             keyboard
%         end
        
    end

end
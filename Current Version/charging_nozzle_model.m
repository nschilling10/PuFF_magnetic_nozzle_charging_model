function [I1, I2, Vfcc, Vcap, R_plasma, V_plasma]=charging_nozzle_model(circInps,plasmaInps,graphDisplay)

% Nathan Schilling
% 05/05/2020
% Takes parameters regarding the plasma and circuit inputs and outputs 

%units
eV = 11605;
keV = 1000*eV;

%constants
mu0 = 4*pi*1e-7;
Ru=8314.5;
e=1.602176634e-19; %[e]=J/eV
k_b=e*8.617333262145e-5; %[k_b]=eV/K

%make sure winding is backward.
% MT = -abs(MT);

% Assign inputs
%primary circuit (1) inputs
L1 = circInps.L1;
R1 = circInps.R1;

%secondary circuit (2) inputs including reactor capacitance C
L2 = circInps.L2;
R2 = circInps.R2;
C = circInps.C;  

%transformer inputs
mu_r = circInps.mu_r;

rT = circInps.rT;
AT = rT^2 * pi;
NT1 = circInps.NT1;
NT2 = circInps.NT2;
lT1 = circInps.lT1;
lT2 = circInps.lT2;

%flux compression generator inputs
Nfcc = circInps.Nfcc;
Rfcc = circInps.Rfcc;

%calculations from inputs
LT1 = mu0 * mu_r * AT*NT1^2/lT1;
LT2 = mu0 * mu_r * AT*NT2^2/lT2;
k=circInps.k;
MT = k*sqrt(LT1*LT2);

%plasma inputs
m_propellant = plasmaInps.m_propellant;
MW = plasmaInps.MW;
Rgas = Ru/MW;
% Temperature must be input in terms of keV
T0 = plasmaInps.T0*keV;
gamma = plasmaInps.g;

%-----------initial conditions
I10 = circInps.I10;
I20 = circInps.I20;
%     V10 = circInps.V10;
V20 = circInps.V20;

%calculation of secondary initial conditions
taufcc = Rfcc/(2*sqrt(2*gamma*Rgas*2*T0));
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

    if graphDisplay
        disp('Plasma bounced at nozzle wall');
    end
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
    
    if graphDisplay
        disp('Capacitor reached peak charging voltage');
    end
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
Vcap = y(:,4);

if graphDisplay
    
    energies=calcCompEnergiesFun(t,I1, I2, Vfcc, Vcap, R_plasma, V_plasma, circInps, plasmaInps);
    Bfield=calcBfield(circInps,R_plasma);
    [eta_KE,eta_circ,eta_sys]=efficiencyCalcFun(energies,plasmaInps)
    graphDisplayFun(t,I1, I2, Vfcc, Vcap, R_plasma, V_plasma,circInps,plasmaInps,energies,Bfield);
end


    function dy = odefuns(t,y)
        %
        %y1 = I1
        %y2 = I2
        %y3 = V1
        %y4 = Vfcc
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
        
        %------Circuit equations-----
        
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
        
    end

end
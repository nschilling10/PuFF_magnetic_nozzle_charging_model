function mass=magneticNozzleDesign(circuitInps,I1_vec,rho,XC)
% Takes the circuit input variables and the vector of I1 currents 
% and outputs the mass of the magnetic nozzle for a nozzle made of 
% material rho and having cross-section XC
Nfcc=circuitInps.Nfcc;
Rfcc=circuitInps.Nfcc;
l=Rfcc;
H=(l/Nfcc);
C=2*pi*Rfcc;
S=sqrt(H^2+C^2);
L=Nfcc*S;
Vol=XC*L;
R1=rho*L/XC;
Heat_pwr=R1*I1_vec.^2;

Rad_pwr=sigma*eps*(0.5*L*sqrt(XC/pi))*(T_ref^4); %assumes metal with specified epsilson,... 
% at reference temperature (this temeprature must be the same temeprature ...
% as the temperature used for the resistivity value), and this equation
% also assumes the wire has a circular cross section. And also only half
% the surface area is able to radiate to free space

Net_heating=Heat_pwr-Rad_pwr;
Net_heating(Net_heating>0)=0;

Q=trapz(t,Net_heating);

mass=Q/(c*(T_ref-4));

mass=Rho*Vol; %Rho here is density

end
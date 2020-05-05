clear all
% Universal gas constant
R=8.134; %[R] = 

%Dimentions of thruster ionization section
rad=0.1; %[rad] = m
l=0.1; %[l] = m
%% Take in reference values from paper
%Plasma Temperature in eV
T=20; 
%Number Density in m^-3
n=1e16*(1e6); 
% Funadmental charge
e=1.6021766208e-19;

%Plasma thermal pressure in Pascals. Found using P=n*e*T
p_plas=n*e*T 

% Applied mangetic field from paper in Tesla
B=(23e3)*1e-4;
% Permeability of vaccum in SI units
mu0=4*pi*1e-7;

% Magnetic pressure P=B^2/(2mu0)
p_mag=B^2/(2*mu0)

% Safety factor for how much higher magnetic pressure should be than plasma
% pressure
r=p_mag/p_plas

% Paramaters of ionization capacitor bank used in paper
C=235.2e-6;
V=18e3;
% Stored Energy of ionization system used in paper
E_cap=0.5*C*V^2;
% Kurt suggested using a 0.25 fudge factor to account for ionization
% losses. This represents the total energy delivered to the plasma
E_plas=0.25*E_cap;
% Another fudge factor I came up with to convert the energy delievered to
% the plasma to the thermal energy the plasma actualy had
conversionFactor=e*T/E_plas
%% Our system parameters
% Parameters to describe the power our system can deliver. These are
% guesses
C=1e-9;
V=1e3;
E_cap=0.5*C*V^2;

% efficiency of capacitor to thermal energy in a plasma varies wildly, 
% but 10% is a reasonable first guess
eta=0.1;

% Find the estimated plasma temperature by converting the energy our system
% has in it to delivered plasma energy, then to plasma thermal energy
T_ev=(eta*E_cap)/(e*3*n*R*pi*rad^2*l)

% Assume our plasma will have 4-order magnitude less number density than
% their plasma. Using this, estimate the total magnetic pressure we will 
% need using the thermal pressure to magnetic pressure conversion 
% coefficient from before, and the estimated plasma thermal pressure (n*e*T_ev)
p_mag=r*n*e*T_ev;

% Get the required magnetic field strength to create the above magnetic
% pressure
B_req=sqrt(2*mu0*p_mag)
B_req=100e-4;
% Pick a current (lower is better)
I=10;
% Find the number of turns. We want this to be about 10
N=B_req/(mu0*I)

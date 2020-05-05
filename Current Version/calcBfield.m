function Bfield = calcBfield(circInps,R_plasma)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
mu0=4*pi*1e-7;

A0=pi*circInps.Rfcc^2;
B0 = mu0 * circInps.Nfcc * circInps.I10 /circInps.Rfcc;  %initial magnetic field, this is an approximation using the radius of the coils as a length scale
Phi0 = B0 * A0;  %assuming this is conserved

% ---------- Force on plasma vs time
%area between coil and plasma
Area = A0 - pi*R_plasma.^2;  %I'm assuming the plasma is a cylinder equal to Rp, it's current radius

%updated magnetic field between plasma
Bfield = Phi0./Area;

%calculate the force and acceleration acting on the plasma
Force = -Bfield.^2 /2/mu0 .* Area;

end


clear all
close all

test.I0=5e+07;
test.L0=400e-6;
test.R1=0;

test.graphDisplay=true;
[E_gain,E_circ] = circuitModelFunv5(test)
Gain=E_gain/E_circ
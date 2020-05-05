clear all
close all

test.I0=1e6;
test.L0=400e-6;
test.R1=0;
test.R2=0;
test.a=1;

test.graphDisplay=true;
[E_gain,E_circ] = circuitModelFunv5_5(test)
Gain=E_gain/E_circ
%% Nathan Schilling
% wrapper script to test ability to simply insert flux compression into
% Kurt's PIT thurster model
clear all
close all

inputs.Lp=1e-6;
inputs.Lc=100e-6;
inputs.Rp=0;
inputs.Rc=0;
inputs.C_load=20e-3;
inputs.L0=360e-6;
inputs.I0=2.5e5;

test.Lp=25e-9;
test.Lc=(6.8e-6);
test.I0=0.5e6;
test.L0=1.4e-6;
test.C_load=5e-3;

[E_gen,E_cap]=kurtCircuitModelFun(inputs)
Gain=E_cap/E_gen
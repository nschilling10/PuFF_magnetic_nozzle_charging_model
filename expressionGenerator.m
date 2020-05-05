clear all
close all
syms L_1 M R_1 dL_dt I1 dI1_dt L_2 R_2 V_cap I2 dI2_dt L12 L_0 l a
%% Generator term is negative
% solns=solve(L_0*dI1_dt + (R_1 - dL_dt)*I1 - M*dI2_dt == 0,...
%     L_2*dI2_dt + R_2*I2 - V_cap - M*dI1_dt == 0,dI1_dt,dI2_dt)
%% transformer as series of windings
% solns=solve(R_1*I1 + l*dI1_dt + L_1*dI1_dt - (1/a)*M*dI2_dt == dL_dt*I1 + L_0*dI1_dt,...
%       a*(L_2*dI2_dt + R_2*I2 + V_cap) == M*dI1_dt + M*dI2_dt*(a-(1/a)),dI1_dt,dI2_dt)
%% Kurt circuit equations?
% solns=solve(dL_dt*I1+L_1*dI1_dt+R_1*I1-L_12*dI2_dt == 0,...
%     L_2*dI2_dt+R_2*I2-V_cap-L_12*dI1_dt == 0,dI1_dt,dI2_dt)
%% Generator term is positive
solns=solve(L_0*dI1_dt + (R_1 + dL_dt)*I1 - M*dI2_dt == 0,...
    L_2*dI2_dt + R_2*I2 + V_cap - M*dI1_dt == 0,dI1_dt,dI2_dt)
%% Russian paper equations
% syms L1 L2 U
% solns=solve(L1*dI1_dt +  dL_dt*I1 + R_1*I1 + L12*dI2_dt == 0,...
%     L2*dI2_dt + R_2*I2 + L12*dI1_dt + U == 0,dI1_dt,dI2_dt);
%% Display results
solns.dI1_dt
solns.dI2_dt
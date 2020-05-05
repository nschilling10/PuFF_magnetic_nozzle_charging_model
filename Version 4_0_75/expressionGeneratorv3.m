% Generates ode45 functions based on circuit equations
clear all
close all
syms V_sqig L_T1 L_T2 L_1 L_2 I_1 R_1 I_2 R_2 M dI1_dt dI2_dt V_cap
eqn1= (L_T1+L_1)*dI1_dt+I_1*R_1-V_sqig-M*dI2_dt == 0;
eqn2= (L_T2+L_2)*dI2_dt+I_2*R_2+V_cap-M*dI1_dt == 0;
solns=solve(eqn1, eqn2,dI1_dt, dI2_dt);
solns.dI1_dt
solns.dI2_dt
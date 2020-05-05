% Generates ode45 functions based on circuit equations
clear all
close all
syms L_c dI1_dt I1 R_1 L_1 M_1 dI2_dt l_1 M_2 I4 dM2
syms L_2 I2 R_2 V_cap l_2
syms Lp dI4_dt dL_p eta_l
eqn1= L_c*dI1_dt+I1*R_1+L_1*dI1_dt-M_1*dI2_dt+l_1*dI1_dt-M_2*dI4_dt-I4*dM2 == 0;
eqn2= L_2*dI2_dt+I2*R_2+V_cap+l_2*dI2_dt-M_1*dI1_dt == 0;
eqn3= Lp*dI4_dt+dL_p*I4-M_2*dI1_dt-I1*dM2+eta_l*I4==0;
solns=solve(eqn1, eqn2, eqn3,dI1_dt, dI2_dt, dI4_dt);
solns.dI1_dt
solns.dI2_dt
solns.dI4_dt
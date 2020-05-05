%% Problem 2 Script
% Nathan Schilling
clear all
close all
% Constant def
hr=3600;
day=24*hr;
yr=365.24*day;
month=yr/12;
AU=1.495e11;
g0=9.80665;
% Input params
R=80*AU;
gamma=0.1;
lambda=0.3;
N=6;
% Part A
alpha=4e-5;
Beta=2;
T_e=(3/2)*((Beta*R*sqrt(alpha/(2*N)))./(1-lambda.^(N/2))).^(2/3);
T_e/(month)
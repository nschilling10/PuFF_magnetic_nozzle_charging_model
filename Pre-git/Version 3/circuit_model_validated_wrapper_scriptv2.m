%% Nathan Schilling
% Transformer circuit model wrapper script
% 02/19/19
clear all
close all
format long

test.graphDisplay=false;
test.I0=1e6;
test.L0=400e-6;
test.R1=1;
test.R2=1;

L1_vec=logspace(-9,-3,1e2);
L2_vec=logspace(-9,-3,1e2);
Gain_mat=ones(length(L2_vec),length(L1_vec));
for i=1:length(L1_vec)
    for j=1:length(L2_vec)
        test.L1=L1_vec(i);
        test.L2=L2_vec(j);
        [E_gain,E_circ] = circuitModelFunValidatedv2_0(test);
        Gain_mat(j,i) = E_gain/E_circ;
    end
end
%% Plotting 
figure(11);
surf(L1_vec*1e9,L2_vec*1e6,Gain_mat)
set(gca,'xscale','log')
set(gca,'yscale','log')
set(gca,'zscale','log')
xlabel('{\boldmath $L_1$ } \textbf{(}{$\mu$H}\textbf{)}','interpreter','latex','fontsize',24)
ylabel('{\boldmath $L_2$ } \textbf{(}{\boldmath$\mu$H}\textbf{)}','interpreter','latex','fontsize',24)
zlabel('\textbf{Gain (}{\boldmath$\frac{\Delta E_{cap}}{\Delta E_{in}}$}\textbf{)}','interpreter','latex','fontsize',24)
        
%% Find max value and display case from it
[rows,colInd_vec]=max(Gain_mat);
[val,rowInd]=max(rows);
colInd=colInd_vec(rowInd);
test.L1=L1_vec(rowInd);
test.L2=L2_vec(colInd);
test.graphDisplay=true;
[E_gain,E_circ] = circuitModelFunValidatedv2_0(test)
E_gain/E_circ
Gain_mat(colInd,rowInd)
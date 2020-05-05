%% Nathan Schilling
% Transformer circuit model wrapper script
% 02/19/19
clear all
close all

L1_vec=logspace(-6,-3,25);
Lratio=logspace(0,2,25);
Gain_mat=ones(length(Lratio),length(L1_vec));
for i=1:length(L1_vec)
    for j=1:length(Lratio)
        test.L1=L1_vec(i);
        test.Lratio=Lratio(j);
        test.graphDisplay=false;
        [E_gain,E_circ] = circuitModelFunv2(test);
        Gain_mat(j,i) = E_gain/E_circ;
        if Gain_mat(j,i) <=0
            Gain_mat(j,i) = 0.01;
        end
    end
end
%% Plotting 
figure(4);
surf(L1_vec*1e6,Lratio,Gain_mat)
set(gca,'xscale','log')
set(gca,'yscale','log')
set(gca,'zscale','log')
xlabel('{\boldmath $L_1$ } \textbf{(}{\boldmath$\mu$H}\textbf{)}','interpreter','latex','fontsize',24)
ylabel('{\boldmath$\frac{L_2}{L_1}$}','interpreter','latex','fontsize',24)
zlabel('\textbf{Gain (}{\boldmath$\frac{\Delta E_{cap}}{\Delta E_{gen}}$}\textbf{)}','interpreter','latex','fontsize',24)
        
%% Find max value and display case from it
[rows,colInd_vec]=max(Gain_mat);
[val,rowInd]=max(rows);
colInd=colInd_vec(rowInd);
test.L1=L1_vec(rowInd);
test.Lratio=Lratio(colInd);
test.graphDisplay=true;
[E_gain,E_circ] = circuitModelFunv2(test)
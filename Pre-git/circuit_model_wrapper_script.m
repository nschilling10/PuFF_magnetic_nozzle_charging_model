%% Nathan Schilling
% Transformer circuit model wrapper script
% 02/19/19
clear all
close all

L1_vec=logspace(-6,-3,1e1);
Lratio=logspace(0,2,1e1);
E_gain_mat=ones(length(Lratio),length(L1_vec));
for i=1:length(L1_vec)
    for j=1:length(Lratio)
        test.L1=L1_vec(i);
        test.Lratio=Lratio(j);
        test.graphDisplay=false;
        E_gain_mat(j,i) = circuitModelFun(test);
        if E_gain_mat(j,i) <=0
            E_gain_mat(j,i) = 0.01;
        end
    end
end
%% Plotting 
figure(4);
surf(L1_vec*1e6,Lratio,E_gain_mat)
set(gca,'xscale','log')
set(gca,'yscale','log')
set(gca,'zscale','log')
xlabel('Generator side of transformer $L_1$ inductance ($\mu$H)','interpreter','latex','fontsize',16)
ylabel('$\frac{L_2}{L_1}$','interpreter','latex','fontsize',16)
zlabel('Energy on capacitor at burnout (J)','interpreter','latex','fontsize',16)
        

%% Find max value and display case from it
[rows,colInd_vec]=max(E_gain_mat);
[val,rowInd]=max(rows);
colInd=colInd_vec(rowInd);
test.L1=L1_vec(rowInd);
test.Lratio=Lratio(colInd);
test.graphDisplay=true;
circuitModelFun(test)


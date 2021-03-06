target_vel=32.67;
precision=2;
target_vel=round(sqrt(plasmaInps.g*(8314.5/plasmaInps.MW)*2*plasmaInps.T0*1e3*11605)*1e-3,precision);
Vpf_v=Vpf_tensor(:);
Inds=find(round(Vpf_v*1e-3,precision) == target_vel);
[i,j,k]=ind2sub([length(Nfcc_vec) length(Rfcc_vec) length(I10_vec)],Inds);
Nfcc_work=Nfcc_vec(i);
Rfcc_work=Rfcc_vec(j);
I10_work=I10_vec(k);
figure(20)
scatter3(Nfcc_work,Rfcc_work,I10_work*1e-6)
set(gca,'xscale','log')
set(gca,'yscale','log')
set(gca,'zscale','log')
xlabel('$N_{fcc}$ (\#)','interpreter','latex','fontsize',24)
ylabel('$R_{fcc}$ (m)','interpreter','latex','fontsize',24)
zlabel('$I_{1_{0}}$ (MA)','interpreter','latex','fontsize',24)
title(strcat('Nozzle designs that resulted in $V_{plasma_{f}}$ = ',num2str(target_vel),'km/s +/-0.0005km/s'),'interpreter','latex','fontsize',26)
mu0=4*pi*1e-7;
B_inital=mu0 * Nfcc_work.*I10_work./Rfcc_work;
figure(19)
plot(B_inital)
avg_B_inital=mean(B_inital)
[Mode_B_initial,mode_index_B_inital]=mode(B_inital)
standard_deviation_B_inital=std(B_inital)
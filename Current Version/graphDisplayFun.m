function graphDisplayFunction(t,I1, I2, Vfcc, Vcap, R_plasma, V_plasma,circInps,plasmaInps,energies,Bfield)

%% Function version of the lines that display graphs related to the charging circuit model
%plot the results
figure(1),
plot(t*1e6,I1/1e6), hold on
plot(t*1e6,I2/1e6), hold on
plot(t*1e6,Vfcc/1e6), hold off
xlabel('t (\mus)');
ylabel('I (MA), V(MV)')
lh = legend('$$I_1$$','$$I_2$$','$$\tilde{V_{fcc}}$$','location','southeast'); legend boxoff;
lh.Interpreter = 'latex';
grid on

figure(2),
plot(t*1e6,Vcap/1e3)
xlabel('t (\mus)');
ylabel('$V_{cap}$ (kV)','interpreter','latex')
grid on
% lh = legend('Recharge Circuit (kV)','Flux Compression (MV)'); legend boxoff;

figure(3),
plot(t*1e6,R_plasma/(circInps.Rfcc/2))
xlabel('t (\mus)');
ylabel('R_{plasma}/R_{coil}')
grid on

figure(4),
plot(t*1e6,V_plasma/1e3)
xlabel('t (\mus)');
ylabel('V_{plasma} (km/s)')
grid on

figure(5),
plot(t*1e6,Bfield)
xlabel('t (\mus)');
ylabel('B (T)')
grid on

%energy in plasma and circuit
figure(6),
semilogy(t*1e6,energies.E_fcc, ...
     t*1e6,energies.E_plasma, ...
     t*1e6,energies.E_cap, ...
     t*1e6,energies.E_therm,...
     t*1e6,energies.E_R1,...
     t*1e6,energies.E_R2,...
     t*1e6,energies.E_tot,...
     t*1e6,energies.E_tot_2)
xlabel('t (\mus)');
ylabel('E (MJ)')
grid on
lh6 = legend('Coil','Plasma','Cap','Thermal','R_1','R_2','Total','Cap_Plas_Coil'); legend boxoff
ylim([1e-3 1e8])
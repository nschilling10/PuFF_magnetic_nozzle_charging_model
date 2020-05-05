%% Nathan Schilling
% Seed current & trade study wrapper script
% 06/28/19
clear all
close all

test.L1=20e-6;
test.L2=1e-6;
test.L0=3.5e-4;
I0_vec=logspace(0,8,1e3);
test.graphDisplay=false;
Gain_vec=zeros(1,length(I0_vec));
Req_circ_e=zeros(1,length(I0_vec));
for i=1:length(Gain_vec)
    test.I0=I0_vec(i);
    [E_gain,E_circ] = circuitModelFunv4(test);
    Req_circ_e(i)=E_circ;
    Gain_vec(i)=E_gain/E_circ;
end
%% Plot results
figure(7)
loglog(I0_vec,Gain_vec)
xlabel('$I_0$ generator side','interpreter','latex','fontsize',24)
ylabel('\textbf{Gain (}{\boldmath$\frac{\Delta E_{cap}}{\Delta E_{gen}}$}\textbf{)}','interpreter','latex','fontsize',24)

figure(8)
loglog(I0_vec,Req_circ_e,I0_vec,10e6*ones(1,length(I0_vec)))
xlabel('$I_0$ generator side','interpreter','latex','fontsize',24)
ylabel('${\Delta E_{gen}}$','interpreter','latex','fontsize',24)
%% Find max
[val,ind]=max(Gain_vec);
test.I0=I0_vec(ind)
val
test.graphDisplay=true;
[E_gain,E_circ] = circuitModelFunv4(test);
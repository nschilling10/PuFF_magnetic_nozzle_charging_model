%% Nathan Schilling
% Russian Transformer circuit model wrapper script

%test.L1=25e-9;
%test.Lratio=25*(6.8e-6)/test.L1;
test.L1=1e-6;
test.Lratio=0.1;
test.graphDisplay=true;
%[E_gain,E_circ] = circuitModelFunv2_5(test)
%[E_gain,E_circ] = circuitModelFunRussian(test)
[E_gain,E_circ] = circuitModelFunv4(test)
% Plot figures from AutomateTPB_sp_bvaries.m
% Plots the variations in the mean of all the variables when b varies
% Fev 09

clear all
%close all

tmax = 500 ;
acc = 200 ;
i = 1 ;

for b = [0.39:0.01:0.5,0.6:0.1:0.9] %[0.39:0.001:0.6,0.605:0.005:1]

    filename = sprintf('CA_b%g_f09c03d02r00001del01.mat',b);
    load(filename)
    
    acc3 = tmax - acc ;
    x1(i) = mean(vegetation(1,acc3:tmax));
    x2(i) = mean(PatchMax);
    x3(i) = mean(StdPatch);
    x4(i) = mean(SkewnessPatch);
    x5(i) = mean(SpCorre);
    
    i = i + 1 ;
    
end
figure
bb = [0.39:0.01:0.5,0.6:0.1:0.9] ;
subplot(3,1,1)
plot(bb,x1,'k',bb,x4,'g',bb,x5,'m')
legend('vege','skew','corr')
subplot(3,1,2)
plot(bb,x3)
title('Std')
subplot(3,1,3)
plot(bb,x2)
title('PatchMax')

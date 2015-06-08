% Plot figures from AutomateSOC_Vasilis_bvaries.m
% Plots the variations in the variables (value of the last simu) when b varies
% Fev 09

clear all
close all

i = 1 ;

for b = [0.39:0.001:0.6,0.605:0.005:1] % [0.39:0.001:0.6,0.605:0.005:1] 
    
    filename = sprintf('Temporal_b%g_f09c03d02r00001del01.mat',b);
    load(filename)
    
    x1(i) = vegetation(end); 
    x3(i) = TempVar;  % follows 100 cells picked at random in the lattice through time
    y3(i) = VARtemp;  % follows the vegetation cover through time
    x4(i) = TempSkew;
    y4(i) = SKtemp;
    x5(i) = TempCorr; 
    y5(i) = CORtemp; 
    
    i = i + 1 ;
    
end

% x1
% x3
% y3
% x4
% y4
% x5
% y5

bb = [0.39:0.001:0.6,0.605:0.005:1] ;% [0.39:0.001:0.6,0.605:0.005:1] ;

figure
plot(bb,x1,'k',bb,x4,'g')
legend('vege','skew')

figure
plot(bb,x1,'k',bb,y4,'b')
legend('vege','skvege')

figure
plot(bb,x1,'k',bb,x5,'m')
legend('vege','corre')

figure
plot(bb,x1,'k',bb,y5,'b')
legend('vege','correvege')

figure
plot(bb,x1,'k',bb,x3/10,'b')
legend('vege','Std/10')

figure
plot(bb,x1,'k',bb,y3*10,'m')
legend('vege','stdvege*10')

% load CA_b0.9_f09c03d02r00001del01.mat
%     
%    count=size(MoyPowerlaw) ;
%     count = count(1) ;
%     i2 = 1 ;
%     for i1 = 1 : count
%         if MoyPowerlaw(i1,2) > 0
%             MoyPowerlawSansZero(i2,2) = MoyPowerlaw(i1,2) ;
%             MoyPowerlawSansZero(i2,1) = MoyPowerlaw(i1,1) ;
%             i2 = i2 + 1 ;
%         end
%     end
%     figure ;
%     plot(MoyPowerlawSansZero(:,1),MoyPowerlawSansZero(:,2),'k.')


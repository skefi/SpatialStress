% Automate asynchrone
% Modèle de Pascual et al, 2002
% Predator preys
% Programme associé aux fonctions Etat suivant, ChooseRandNeig, NbPreyNeig
% After the transient, saves the patch size distri every 40 steps
% Fev 2009

clear all 
close all

global n z deltat beta1 beta2 s v

tic

z = 4 ;
v = 1 ;
%beta1 = 0.333333 ; % birth rate of preys
beta1 = 0.05 ;
beta2 = 0.1 ; % repro rate of predators which ate a prey
s = 0.15; %0.333333 ; death rate of predators which did not eat a prey

deltat = 1 ; %Pas de temps de l'automate cellulaire
n = 100 ;
tfin = 1000 ; % 10000 ;
acc = 700 ; %7000 ; nb d'iterations sur lesquelles on calcule la distri de taille
%de patchs 
corre = 40 ; % number of steps to skip because there are temporal auto-
% correlation below that  
incretaille = 5 ; % taille d'une classe de taille
itet = 1 ;

%Initialisation de la matrice
%==============================================================
M = rand (n,n) ;
po = 0.5 ; %Proportion of prey
pd = 0.5 ; %Proportion of empty sites; 1-pd-po = proportion of pred
DensitePrey = 0 ; % prey
DensitePred = 0 ; % empty

for i = 1 : n
    for k = 1 : n
        if M(i,k)<po
            M(i,k) = 1 ; %prey
            DensitePrey = DensitePrey + 1 ;
        elseif M(i,k) > pd + po
            M(i,k) = 2 ; %pred
            DensitePred = DensitePred + 1 ;
        else M(i,k) = 3 ; %empty
        end
    end
end
Mo = M ; %Matrice de l'état initial

Mana = (M==1) ;
[SortieMana,NbPatch] = bwlabel(Mana,4);
PatchData = regionprops(SortieMana,'basic');
AllPatch1 = [PatchData.Area];

prey(1) = DensitePrey/(n*n) ;
pred(1) = DensitePred/(n*n) ;


%Progression du système
%==============================================================
tmax = tfin/deltat ;
for t = 0 : tmax

%     black = [0 0 0] ;
%     gray = [0.7 0.7 0.7] ;
%     white = [1 1 1] ;
%     Satecolormap = [black; gray; white] ; % 1=prey=noir 2=pred=grey et 3=empty=white
%     colormap(Satecolormap) ;
%     image(M) ; axis off ;
%     drawnow ;

    [M] = EtatSuivant(M) ;

    %Obtention des densités globales
    DensitePrey = sum(sum((M==1))) ;
    DensitePred = sum(sum((M==2))) ;
    prey(t+2) = DensitePrey/(n*n) ;
    pred(t+2) = DensitePred/(n*n) ;

    %======================================================================
    %Distribution de taille de patch sur les acc derniers pas de temp
    %======================================================================

    if t > tmax-acc
        
        if mod(t,corre) == 0

            Mana = (M==1) ;
            [SortieMana,NbPatch] = bwlabel(Mana,4);
            PatchData = regionprops(SortieMana,'basic');
            AllPatch = [PatchData.Area];

            nbpoints = n*n/incretaille ;
            taille = 0 ;
            s1 = size(AllPatch) ;
            s1 = s1(2) ;
            increi = 1 ;
            xreg = 0 ;

            for i = 1:nbpoints + 1

                NbPatchSuptaille = 0 ;

                for k = 1:s1
                    if ( AllPatch(1,k)>taille & AllPatch(1,k)<=taille+incretaille )
                        NbPatchSuptaille = NbPatchSuptaille + 1 ;
                    end
                end
                xreg(increi,1) = taille + incretaille/2 ; % Size
                xreg(increi,2) = NbPatchSuptaille ; % nombre de patchs qui ont
                %une taille superieure a taille
                taille = taille + incretaille ;
                increi = increi + 1 ;

            end

            xreg = sortrows(xreg) ;
            taillexreg = size(xreg) ;
            lixreg = taillexreg(1) ;
            colxreg = taillexreg(2) ;
            for i1 = 1 : lixreg
                for i2 = 1 : colxreg
                    Powerlaw(i1,i2,itet) = xreg(i1,i2) ;
                end
            end
            
            % Suivi du patch le plus grand
            PatchMax(1,itet) = max(AllPatch) ;

            itet = itet+1 ;
        end
    end % if t>

end %Fin boucle sur t

%==========================================================================
%Tracé des courbes
%==========================================================================
toc
acc2 = itet-1 ;
MoyPowerlaw = Powerlaw(:,1,1) ;
MoyPowerlaw(:,2) = 0 ;
for i = 1 : acc2
    MoyPowerlaw(:,2) = MoyPowerlaw(:,2) + Powerlaw(:,2,i) ;
end
MoyPowerlaw(:,2) = MoyPowerlaw(:,2)./acc2 ;

count=size(MoyPowerlaw) ;
count = count(1) ;
i2 = 1 ;
for i1 = 1 : count
    if MoyPowerlaw(i1,2) > 0
        MoyPowerlawSansZero(i2,2) = MoyPowerlaw(i1,2) ;
        MoyPowerlawSansZero(i2,1) = MoyPowerlaw(i1,1) ;
        i2 = i2 + 1 ;
    end
end
figure ;
plot(MoyPowerlawSansZero(:,1),MoyPowerlawSansZero(:,2),'k.')

figure ;
black = [0 0 0] ;
gray = [0.7 0.7 0.7] ;
white = [1 1 1] ;
Satecolormap = [black; white; gray] ; % 1=noir 2=gris et 3=white
colormap(Satecolormap) ;
image(M) ; axis off ;

acc3 = tmax - acc ;
mean(prey(1,acc3:tmax))
mean(pred(1,acc3:tmax))
mean(PatchMax)

figure ;
plot(prey)
hold on
plot(pred,'k')

% save Pasc_Skip40_b10333b201s08.mat MoyPowerlaw M prey pred
% MoyPowerlawSansZero(:,1)
 

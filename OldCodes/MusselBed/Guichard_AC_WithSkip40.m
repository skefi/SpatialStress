% Automate asynchrone
% Modèle de Guichard et al, 2003
% Mussel beds
% Programme associé aux fonctions EtatSuivant1, MatriceVoisins2 
% After the transient, saves the patch size distri every 40 steps
% Dec 2008

% PatchDistri2=PatchDistri';
% dlmwrite('r015_rough.txt',PatchDistri2,'delimiter','\t')

clear all 
close all

global n z deltat d alpha2 delta0 alpha0

tic

z = 4 ;
d = 1 ; %
alpha2 = 0.8 ; % productivity
delta0 = 0.03 ; % basic rate from 2 to 0
alpha0 = 0.1 ; % wave exposure

deltat = 1 ; %Pas de temps de l'automate cellulaire
n = 50   %100 ;
tfin = 100%10000
acc = 70%7000 ; %nb d'iterations sur lesquelles on calcule la distri de taille
%de patchs 
corre = 4%40 ; % number of steps to skip because there are temporal auto-
% correlation below that  
incretaille = 5 ; % taille d'une classe de taille
itet = 1 ;

%Initialisation de la matrice
%==============================================================
M = rand (n,n) ;
po = 0.2 ; %Proportion of mussels
pd = 0.6 ; %Proportion of disturbed sites; 1-pd-po = proportion de sites vides
DensiteOcc = 0 ;
DensiteDe = 0 ;

for i=1:n
    for k = 1:n
        if M(i,k)<po
            M(i,k) = 1 ; %mussel
            DensiteOcc = DensiteOcc + 1 ;
        elseif M(i,k)>pd + po
            M(i,k) = 2 ; %empty
        else M(i,k) = 3 ; %disturbed
            DensiteDe = DensiteDe + 1 ;
        end
    end
end
Mo = M ; %Matrice de l'état initial

Mana = (M==1) ;
[SortieMana,NbPatch] = bwlabel(Mana,4);
PatchData = regionprops(SortieMana,'basic');
AllPatch1 = [PatchData.Area];

mussel(1) = DensiteOcc/(n*n) ;
disturbed(1) = DensiteDe/(n*n) ;

%Initialisation de la matrice V qui à chaque individu associe son
%nombre de voisins occupés
V2 = MatriceVoisinsMoules(M) ;
V0 = MatriceVoisinsDisturbed(M) ;

% %Initiation des densités locales
% qplus(1) = 0 ; %q+|+
% qmoins(1) = 0 ; %q+|-
% for i = 1:n
%     for k = 1:n
%         if M(i,k) == 1 % mussel
%             qplus(1) = qplus(1) + V2(i,k)/4 ;
%         elseif M(i,k) == 3 %disturbed
%             qmoins(1) = qmoins(1) + V2(i,k)/4 ;
%         end
%     end
% end
% qplus(1) = qplus(1)/DensiteOcc ; %q+|+
% qmoins(1) = qmoins(1)/DensiteDe ; %q+|-
% qmoins(1) = qmoins(1)*DensiteDe/DensiteOcc ; %q-|+

%Progression du système
%==============================================================
tmax = tfin/deltat ;
for t = 0 : tmax

%     black = [0 0 0] ;
%     gray = [0.7 0.7 0.7] ;
%     white = [1 1 1] ;
%     Satecolormap = [black; white; gray] ; % 1=mussel=noir 2=empty=white et 3=perturbed=grey
%     colormap(Satecolormap) ;
%     image(M) ; axis off ;
%     drawnow ;

    [M,V0,V2,DensiteOcc,DensiteDe] = EtatSuivant(M,V0,V2,DensiteOcc,DensiteDe) ;

    %Obtention des densités globales
    mussel(t+2) = DensiteOcc/(n*n) ;
    disturbed(t+2) = DensiteDe/(n*n) ;

%     %Obtention des densités locales
%     qplus(t+2) = 0 ; %q+|+
%     qmoins(t+2) = 0 ; %q+|-
%     for i = 1:n
%         for k = 1:n
%             if M(i,k) == 1 %Végétation
%                 qplus(t+2) = qplus(t+2) + V2(i,k)/4 ;
%             elseif M(i,k) ==3 %Dégradé
%                 qmoins(t+2) = qmoins(t+2) + V2(i,k)/4 ;
%             end
%         end
%     end

%     if DensiteOcc>0.00000001
%         qplus(t+2) = qplus(t+2)/DensiteOcc ; %q+|+
%         qmoins(t+2) = qmoins(t+2)*DensiteDe/DensiteOcc ; %q-|+
%     else
%         qplus(t+2) = 0 ;
%         qmoins(t+2) = 0 ;
%     end
%     if DensiteDe>0.00000001
%         qmoins(t+2) = qmoins(t+2)/DensiteDe ; %q+|-
%     else
%         qmoins(t+2) = 0 ;
%     end

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

figure
plot(mussel)

acc3 = tmax - acc ;
mean(mussel(1,acc3:tmax))
mean(PatchMax)

% figure ;
% plot(mussel)
% hold on
% plot(disturbed,'k')

res(:,1) = MoyPowerlawSansZero(:,1) ; %size
res(:,2) = MoyPowerlawSansZero(:,2) ; %number
res(:,3) = log(MoyPowerlawSansZero(:,1)) ; %log size
res(:,4) = log(MoyPowerlawSansZero(:,2)) ; %log number

%save Gui_Skip40_do0041a208ao01.mat MoyPowerlaw M mussel
%dlmwrite('Gui_Skip40_do00410a208ao01.txt',res,'delimiter','\t')
 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


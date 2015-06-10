% Automate asynchrone, TPB
% Project with Vasilis
% Programme associé aux fonctions EtatSuivant1, MatriceVoisins2
% Fev 09

clear all 
close all

global n m r f b d j del z deltat

tic

% Para values of the Nature paper
% -------------------------------------------------------------------------
% for m = 0.1, discontinuous transition
% extinction of the system for b<0.39
z = 4 ;
f = 0.3 ;
d = 0.2 ;
j = 0.0001 ;
m = 0.1 ;
r = 0.9 ;
del = 0.1 ;
deltat = 0.1 ; %Pas de temps de l'automate cellulaire
n = 100 ;
tfin = 1000 ; % tot nb of time steps = 50/deltat = 500
acc = 200 ; %nb d'iterations sur lesquelles on calcule la distri de taille
%de patchs
incretaille = 5 ; % taille d'une classe de taille

%for b = [0.39:0.001:0.6,0.605:0.005:1]
    
    b = 0.8
    %Initialisation de la matrice
    %==============================================================
    M = rand (n,n) ;
    po = 0.6 ; %Proportion de végétation
    pd = 0.2 ; %1-pd-po = proportion de sites vides
    DensiteOcc = 0 ;
    DensiteDe = 0 ;

    for i=1:n
        for k = 1:n
            if M(i,k)<po
                M(i,k) = 1 ; %végétation
                DensiteOcc = DensiteOcc + 1 ;
            elseif M(i,k)>pd + po
                M(i,k) = 2 ; %vide
            else M(i,k) = 3 ; %désert
                DensiteDe = DensiteDe + 1 ;
            end
        end
    end

    vegetation(1) = DensiteOcc/(n*n) ;
    degrade(1) = DensiteDe/(n*n) ;

    %Initialisation de la matrice V qui à chaque individu associe son
    %nombre de voisins occupés
    V = MatriceVoisins2(M) ;

    %Initiation des densités locales
    qplus(1) = 0 ; %q+|+
    qmoins(1) = 0 ; %q+|-
    for i = 1:n
        for k = 1:n
            if M(i,k) == 1 %Végétation
                qplus(1) = qplus(1) + V(i,k)/4 ;
            elseif M(i,k) ==3 %Dégradé
                qmoins(1) = qmoins(1) + V(i,k)/4 ;
            end
        end
    end
    qplus(1) = qplus(1)/DensiteOcc ; %q+|+
    qmoins(1) = qmoins(1)/DensiteDe ; %q+|-
    qmoins(1) = qmoins(1)*DensiteDe/DensiteOcc ; %q-|+

    %Progression du système
    %==============================================================
    tmax = tfin/deltat ;
    itet = 1 ;
    for t = 0 : tmax

        [M,V,DensiteOcc,DensiteDe] = EtatSuivant1(M,V,DensiteOcc,DensiteDe) ;

        %Obtention des densités globales
        vegetation(t+2) = DensiteOcc/(n*n) ;
        degrade(t+2) = DensiteDe/(n*n) ;

        %Obtention des densités locales
        qplus(t+2) = 0 ; %q+|+
        qmoins(t+2) = 0 ; %q+|-
        for i = 1:n
            for k = 1:n
                if M(i,k) == 1 %Végétation
                    qplus(t+2) = qplus(t+2) + V(i,k)/4 ;
                elseif M(i,k) ==3 %Dégradé
                    qmoins(t+2) = qmoins(t+2) + V(i,k)/4 ;
                end
            end
        end

        if DensiteOcc>0.00000001
            qplus(t+2) = qplus(t+2)/DensiteOcc ; %q+|+
            qmoins(t+2) = qmoins(t+2)*DensiteDe/DensiteOcc ; %q-|+
        else
            qplus(t+2) = 0 ;
            qmoins(t+2) = 0 ;
        end
        if DensiteDe>0.00000001
            qmoins(t+2) = qmoins(t+2)/DensiteDe ; %q+|-
        else
            qmoins(t+2) = 0 ;
        end

        %======================================================================
        %Distribution de taille de patch sur les acc derniers pas de temp
        %======================================================================

        if t > tmax-acc
            for i=1:n
                for k=1:n
                    if M(i,k) == 1
                        Mana(i,k) = 1 ;
                    else
                        Mana(i,k) = 0 ;
                    end
                end
            end
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
            %VariancePatch(1,itet) = var(AllPatch) ;
            StdPatch(1,itet) = std(AllPatch) ;
            SkewnessPatch(1,itet) = skewness(AllPatch) ;
            SpCorre(1,itet) = qplus(t+2)/vegetation(t+2) ; % spatial correlation from TPB
            
            itet = itet+1 ;

        end % if t>

    end %Fin boucle sur t

    %==========================================================================
    %Tracé des courbes
    %==========================================================================
    toc
    MoyPowerlaw = Powerlaw(:,1,1) ;
    MoyPowerlaw(:,2) = 0 ;
    for i = 1 : itet-1
        MoyPowerlaw(:,2) = MoyPowerlaw(:,2) + Powerlaw(:,2,i) ;
    end
    MoyPowerlaw(:,2) = MoyPowerlaw(:,2)./acc ;    
    
    %filename = sprintf('CA_b%g_f09c03d02r00001del01.mat',b)
    %save(filename,'MoyPowerlaw', 'M', 'vegetation', 'PatchMax', 'StdPatch', 'SkewnessPatch', 'SpCorre')
    
    hist(AllPatch,20);
    %save AllPatchb08.mat AllPatch
    
    %clear AllPatch PatchData xreg Powerlaw PatchMax StdPatch SkewnessPatch SpCorre MoyPowerlaw

%end


%     count=size(MoyPowerlaw) ;
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

%     figure ;
%     black = [0 0 0] ;
%     gray = [0.7 0.7 0.7] ;
%     white = [1 1 1] ;
%     Satecolormap = [black; gray; white] ; % 1=noir 2=gris et 3=white
%     colormap(Satecolormap) ;
%     image(M) ; axis off ;
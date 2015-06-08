% Automate asynchrone, TPB
% Programme associé aux fonctions EtatSuivant1, MatriceVoisins2
% Calculate temporal indicators
% March 09

% ax1 = gca;
% set(ax1);
% ax2 = axes('Position',get(ax1,'Position'),'YAxisLocation','right','Color','none','YColor',[0.5 0.5 0.5]);
% X(:,1) = [0.2; 0.3; 0.4; 0.42; 0.44; 0.46; 0.48; 0.5; 0.52; 0.54; 0.56; 0.58; 0.6; 0.65; 0.7; 0.8; 0.9; 1] ;
% X(:,2) = [0; 0; 76.3895; 105.52; 161.191; 223.7945; 329.3915; 595.236; 974.8; 1682.8; 3270.4; 4470.7; 5448.8; 6253.4; 6827.4; 7461.2; 7891.1; 8179.6] ;
% % plot(X1(:,1),X1(:,2)/7000)
% hl1 = line(X(:,1),X(:,2)/10000, 'Color',[0.5 0.5 0.5], 'Linestyle','--', 'LineWidth', 3, 'Parent', ax2);

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

for b = [0.593:0.001:0.6,0.605:0.005:1]

    b
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
    AllPatch1 = [PatchData.Area];

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

    % Pick 100 random cells in the lattice
    ncells = 100 ;
    randomCells = ceil(n*n*rand(sqrt(ncells))) ; % gives a matrix of ncells*ncells
    % cells that each contain an integer between 1 and n*n (cells picked in the
    % lattice)
    SubM=M(randomCells) ;
    MM=SubM(:);
    TimeData = MM ;

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
            VariancePatch(1,itet) = var(AllPatch) ;
            SkewnessPatch(1,itet) = skewness(AllPatch) ;
            SpCorre(1,itet) = qplus(t+2)/vegetation(t+2) ; % spatial correlation from TPB

            % pick and follow 100 random cells
            SubM = M(randomCells) ; % extracts the matrix of only the random cells
            MM = SubM(:) ; % transforms the matrix into a vector
            TimeData(1:ncells,itet) = MM ; % saves the vector in TimeData

            itet = itet+1 ;

        end % if t>

    end %Fin boucle sur t

    %==========================================================================
    % Temporal indicators
    %==========================================================================

    % Temporal statistics on the vegetation cover
    MEANtemp = mean(vegetation(tmax-acc+1:tmax)); % calculate mean of each timeseries
    VARtemp = std(vegetation(tmax-acc+1:tmax)); % calculate SD of each timeseries
    SKtemp = skewness(vegetation(tmax-acc+1:tmax)); % calculate skewness of each timeseries
    [w, CORtemp] = arfit(vegetation(tmax-acc+1:tmax)', 1, 1); % w: weight/error

    % Temporal stat on 100 randomly picked cells
    %size(TimeData)
    TimeData = (TimeData == 1) ;
    ProbaSucces = sum(TimeData,2)/acc ;
    ProbaEchec = max(0,1 - ProbaSucces) ;
    VarianceData = ProbaSucces .* acc .*ProbaEchec ;

    for count = 1:ncells
        if VarianceData(count) == 0
            SkewnessData(count,1) = 0 ;
        else
            SkewnessData(count,1) = (1-2*ProbaSucces(count))/(sqrt(VarianceData(count))) ;
        end
    end

    corr = zeros(ncells,1) ;
    for line = 1 : ncells
        for col = 2 : acc
            if (TimeData(line,col) == 1 & TimeData(line,col-1) == 1)
                corr(line) = corr(line) + 1 ;
            end
        end
        tt = sum(TimeData(line,:)) ;
        if tt == 0
            corr(line) = 0 ;
        else
            corr(line) = corr(line)/tt ;  % should it be tt - 1 if there is more that 2 occupied sites?
        end
    end

    TempVar = mean(VarianceData) ;
    TempSkew = mean(SkewnessData) ;
    TempCorr = mean(corr) ;

    %==========================================================================
    %Tracé des courbes
    %==========================================================================

    MoyPowerlaw = Powerlaw(:,1,1) ;
    MoyPowerlaw(:,2) = 0 ;
    for i = 1 : itet-1
        MoyPowerlaw(:,2) = MoyPowerlaw(:,2) + Powerlaw(:,2,i) ;
    end
    MoyPowerlaw(:,2) = MoyPowerlaw(:,2)./acc ;

    filename = sprintf('Temporal_b%g_f09c03d02r00001del01.mat',b) ;
    save(filename,'MoyPowerlaw', 'M', 'vegetation','TempVar','TempSkew','TempCorr','VARtemp','SKtemp','CORtemp')

    clear AllPatch PatchData xreg Powerlaw PatchMax StdPatch SkewnessPatch SpCorre MoyPowerlaw ProbaSucces VarianceData SkewnessData Corr

end
toc
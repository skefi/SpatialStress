% Predator prey
% Fonction permettant le calcul de l'état système à t+1 sachant t
% Automate synchrone
% Associé au programme principal: Pascual_AC_WithSkip40.m
% Associé aux fonctions: ChooseRandNeig, NbPreyNeig

function[N] = EtatSuivant(M)

global n z deltat beta1 beta2 s v

N = M ;

for i = 1 : n*n

    %Identification de la cellule concernée
    %----------------------------------------------------------------------
    S = floor(rand*(n*n-1)+1) ; %N° de la cellule, nb entre 1 et n*n
    % S = n * (Si-1) + Sj

    Sj = mod(S,n) ; %Colonne correspondante
    if Sj == 0
        Sj = n ;
        Si = floor(S/n) ; %Ligne correspondante
    else
        Si = floor(S/n) + 1; %Ligne correspondante
    end

    if Si == 0
        Si = n ;
    end

    %Modification de l'état?
    %----------------------------------------------------------------------
    test = rand ;

    if N(Si,Sj) == 1 %Occupe par une proie

        x = ChooseRandNeig(Si,Sj) ; % function which gives the coordinates
        % x = [i,j] of one of the 4 neighbors randomly chosen
        x1 = x(1) ;
        x2 = x(2) ;
        if M(x1,x2) == 3 % the chosen neighbor is empty
            if test < beta1 * deltat
                N(x1,x2) = 1 ; % the chosen neighbor becomes occupied by a prey (repro)
                %DensitePrey = DensitePrey + 1 ;
                %DensiteEmpty = DensiteEmpty - 1 ;
            end
        end

    elseif N(Si,Sj) == 2 %Occupe par un predateur
        nprey = NbPreyNeig(M,Si,Sj) ; % coordinates of the prey neighbors
        TotPrey = sum(nprey(1,:)) ; % tot nb of prey neignbors

        if TotPrey > 0 % the pred has at least one prey neighbor

            select = floor(rand*TotPrey) ; % the pred chooses one prey at random
            if select == 0
                select = TotPrey ;
            end
            Si2 = nprey(2,select) ; %prey site replaced by pred
            Sj2 = nprey(3,select) ;
            N(Si2,Sj2) = 2 ;
            %DensitePrey = DensitePrey - 1 ;

            if test < (1-beta2) * deltat % pred does not reproduce
                N(Si,Sj) = 3 ; % the original site of the pred becomes empty
                %DensiteEmpty = DensiteEmpty + 1 ;
            end

        else % no prey in the neighborhood

            if test < s
                N(Si,Sj) = 3 ; % pred starves
                %DensiteEmpty = DensiteEmpty + 1 ;
            end
            
        end % end cases for pred
    
    end % end boucle on states
    
end % end boucle for

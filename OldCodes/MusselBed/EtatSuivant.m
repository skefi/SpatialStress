% Fonction permettant le calcul de l'état système à t+1 sachant t
% Automate synchrone

function [N,W0,W2,DensiteOcc2,DensiteDe2] = EtatSuivant(M,V0,V2,DensiteOcc,DensiteDe)

N = M ;
W0 = V0 ;
W2 = V2 ;
DensiteDe2 = DensiteDe ;
DensiteOcc2 = DensiteOcc ;

global n z deltat d alpha2 delta0 alpha0

for i = 1 : n*n

    %Identification de la cellule concernée
    S = floor(rand*n*n) ; %N° de la cellule
    Si = floor(S/n) ; %Ligne correspondante
    Sj = mod(S,n) ; %Colonne correspondante
    if Si == 0
        Si = n ;
    end
    if Sj == 0
        Sj = n ;
    end

    %Modification de l'état?
    test = rand ;
    if N(Si,Sj) == 1 %Occupe par une moule

        if W0(Si,Sj) > 0 % if the site has at least one disturbed neighbour
            mu0 = 1;
        else
            mu0 = 0;
        end
        if (delta0 +  mu0*alpha0)*deltat > 1
            disp('Transition mussel --> disturbed certaine. Diminuer deltat')
        end
        if test < (delta0 +  mu0*alpha0)*deltat
            N(Si,Sj) = 3 ; % mussel site becomes a disturbed site

            % Update mussel and disturbed matrices
            z = mod(Sj-1+n,n) ;
            if z == 0
                z = n ;
            end
            W2(Si,z) = W2(Si,z) - 1 ;
            W0(Si,z) = W0(Si,z) + 1 ;
            z = mod(Si-1+n,n) ;
            if z == 0
                z = n ;
            end
            W2(z,Sj) = W2(z,Sj) - 1 ;
            W0(z,Sj) = W0(z,Sj) + 1 ;
            z = mod(Sj+1,n) ;
            if z == 0
                z = n ;
            end
            W2(Si,z) = W2(Si,z) - 1 ;
            W0(Si,z) = W0(Si,z) + 1 ;
            z = mod(Si+1,n) ;
            if z == 0
                z = n ;
            end
            W2(z,Sj) = W2(z,Sj) - 1 ;
            W0(z,Sj) = W0(z,Sj) + 1 ;

            DensiteOcc2 = DensiteOcc2 - 1 ;
            DensiteDe2 = DensiteDe2 + 1 ;
        end

    elseif N(Si,Sj) == 3 %Disturbed site

        % if test < d*deltat
        N(Si,Sj) = 2 ; %Regeneration

        % Update disturbed matrix
        z = mod(Sj-1+n,n) ;
        if z == 0
            z = n ;
        end
        W0(Si,z) = W0(Si,z) - 1 ;
        z = mod(Si-1+n,n) ;
        if z == 0
            z = n ;
        end
        W0(z,Sj) = W0(z,Sj) - 1 ;
        z = mod(Sj+1,n) ;
        if z == 0
            z = n ;
        end
        W0(Si,z) = W0(Si,z) - 1 ;
        z = mod(Si+1,n) ;
        if z == 0
            z = n ;
        end
        W0(z,Sj) = W0(z,Sj) - 1 ;

        DensiteDe2 = DensiteDe2 - 1 ;
        %end

    elseif test < ( alpha2*W2(Si,Sj)/z )*deltat % the site is empty (2)

        N(Si,Sj) = 1 ; %Recolonisation

        % Update matrix mussel
        z = mod(Sj-1+n,n) ;
        if z == 0
            z = n ;
        end
        W2(Si,z) = W2(Si,z) + 1 ;
        z = mod(Si-1+n,n) ;
        if z == 0
            z = n ;
        end
        W2(z,Sj) = W2(z,Sj) + 1 ;
        z = mod(Sj+1,n) ;
        if z == 0
            z = n ;
        end
        W2(Si,z) = W2(Si,z) + 1 ;
        z = mod(Si+1,n) ;
        if z == 0
            z = n ;
        end
        W2(z,Sj) = W2(z,Sj) + 1 ;

        DensiteOcc2 = DensiteOcc2 + 1 ;

    end %Fin de la boucle d'état de la cellule (1, 2 ou 3)
end %Fin de la boucle for (n*n tirages de cellules)



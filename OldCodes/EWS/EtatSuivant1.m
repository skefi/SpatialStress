% Fonction permettant le calcul de l'état système à t+1 sachant t
% Automate synchrone

function [N,W,DensiteOcc2,DensiteDe2] = EtatSuivant (M,V,DensiteOcc,DensiteDe)

N = M ;
W = V ;
DensiteDe2 = DensiteDe ;
DensiteOcc2 = DensiteOcc ;

global n m r f b d j del deltat

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
    if N(Si,Sj) == 1 %Occupe par de la végétation
        if test < m*deltat
            N(Si,Sj) = 2 ; %Mortalite vegetation
            z = mod(Sj-1+n,n) ;
            if z == 0
                z = n ;
            end
            W(Si,z) = W(Si,z) - 1 ;
            z = mod(Si-1+n,n) ;
            if z == 0
                z = n ;
            end
            W(z,Sj) = W(z,Sj) - 1 ;
            z = mod(Sj+1,n) ;
            if z == 0
                z = n ;
            end
            W(Si,z) = W(Si,z) - 1 ;
            z = mod(Si+1,n) ;
            if z == 0
                z = n ;
            end
            W(z,Sj) = W(z,Sj) - 1 ;
            DensiteOcc2 = DensiteOcc2 - 1 ;
        end
        
    elseif N(Si,Sj) == 3 %Occupe par un site dégradé
        if (j + r*W(Si,Sj)/4)*deltat > 1
            disp('Transition - --> o certaine. Diminuer deltat')
        end
        if test < (j + r*W(Si,Sj)/4)*deltat 
            N(Si,Sj) = 2 ; %Regeneration 
            DensiteDe2 = DensiteDe2 - 1 ;
        end
   
    elseif test < ( (del*DensiteOcc2/(n*n)+(1-del)*W(Si,Sj)/4)*(b-f*DensiteOcc2/(n*n)) )*deltat 
         if ( (del*DensiteOcc2/(n*n)+(1-del)*W(Si,Sj)/4)*(b-f*DensiteOcc2/(n*n))+ d )*deltat > 1
             disp('Transition o --> o impossible. Diminuer deltat')
         end
         N(Si,Sj) = 1 ; %Recolonisation
         z = mod(Sj-1+n,n) ;
         if z == 0
             z = n ;
         end
         W(Si,z) = W(Si,z) + 1 ;
         z = mod(Si-1+n,n) ;
         if z == 0
             z = n ;
         end
         W(z,Sj) = W(z,Sj) + 1 ;
         z = mod(Sj+1,n) ;
         if z == 0
             z = n ;
         end
         W(Si,z) = W(Si,z) + 1 ;
         z = mod(Si+1,n) ;
         if z == 0
             z = n ;
         end
         W(z,Sj) = W(z,Sj) + 1 ;
         DensiteOcc2 = DensiteOcc2 + 1 ;
          
     elseif  test > 1 - d*deltat 
         N(Si,Sj) = 3 ; %Degradation
         DensiteDe2 = DensiteDe2 + 1 ;
     end %Fin de la boucle d'état de la cellule (1, 2 ou 3)
 end %Fin de la boucle for (n*n tirages de cellules)
 


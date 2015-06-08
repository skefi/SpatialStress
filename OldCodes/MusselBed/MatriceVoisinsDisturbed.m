% Fonction qui construit la matrice V qui à chaque site associe son
% Nombre de voisins disturbed
% Automate synchrone

function[V0] = MatriceVoisinsDisturbed(M)

global n

V0 = zeros(n,n) ;

for i = 1 : n
    for j = 1 : n
        z = mod(j-1+n,n) ;
        if z==0
            z=n;
        end
        if M(i,z) == 3
            V0(i,j) = V0(i,j)+1 ;
        end
        z = mod(i-1+n,n) ;
        if z==0
            z=n;
        end
        if M(z,j) == 3
            V0(i,j) = V0(i,j)+1 ;
        end
        z = mod(j+1,n) ;
        if z==0
            z=n;
        end
        if M(i,z) == 3
            V0(i,j) = V0(i,j)+1 ;
        end
        z = mod(i+1,n);
        if z==0
            z=n;
        end
        if M(z,j) == 3
            V0(i,j) = V0(i,j)+1 ;
        end
    end
end


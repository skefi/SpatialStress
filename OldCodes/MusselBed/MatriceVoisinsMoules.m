% Fonction qui construit la matrice V qui à chaque site associe son
% Nombre de voisins occupés par des moules
% Automate synchrone

function[V2] = MatriceVoisinsMoules(M)

global n

V2 = zeros(n,n) ;

for i = 1 : n
    for j = 1 : n
        z = mod(j-1+n,n) ;
        if z==0
            z=n;
        end
        if M(i,z) == 1
            V2(i,j) = V2(i,j)+1 ;
        end
        z = mod(i-1+n,n) ;
        if z==0
            z=n;
        end
        if M(z,j) == 1
            V2(i,j) = V2(i,j)+1 ;
        end
        z = mod(j+1,n) ;
        if z==0
            z=n;
        end
        if M(i,z) == 1
            V2(i,j) = V2(i,j)+1 ;
        end
        z = mod(i+1,n);
        if z==0
            z=n;
        end
        if M(z,j) == 1
            V2(i,j) = V2(i,j)+1 ;
        end
    end
end


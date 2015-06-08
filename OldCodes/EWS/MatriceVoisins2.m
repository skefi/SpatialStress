%Fonction qui construit la matrice V qui à chaque site associe son
% Nombre de voisins occupés
% Automate synchrone

function[V] = MatriceVoisins(M)

global n

V = zeros(n,n) ;

%Centre
for i = 1 : n
    for j = 1 : n
        z = mod(j-1+n,n) ;
        if z==0
            z=n;
        end
        if M(i,z) == 1
            V(i,j) = V(i,j)+1 ;
        end
        z = mod(i-1+n,n) ;
        if z==0
            z=n;
        end
        if M(z,j) == 1
            V(i,j) = V(i,j)+1 ;
        end
        z = mod(j+1,n) ;
        if z==0
            z=n;
        end
        if M(i,z) == 1
            V(i,j) = V(i,j)+1 ;
        end
        z = mod(i+1,n);
        if z==0
            z=n;
        end
        if M(z,j) == 1
            V(i,j) = V(i,j)+1 ;
        end
    end
end


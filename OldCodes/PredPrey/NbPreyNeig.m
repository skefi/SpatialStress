% Fonction qui associe à un site donné le nombre de voisins occupés par des
% proies et leur coordonnees dans le reseau
% Associated with the pg Pascual_AC_WithSkip40.m

% The total nb of prey neighbors is the sum of the first line
% sum(nbprey(1,:))
% Chaque colonne correspond aux coordonnées d'une proie (i:ligne 2, j:ligne 3) 

function[nprey] = NbPreyNeig(M,Si,Sj)

global n

% M = [ 1 2 3 ; 2 1 2 ; 3 3 3] 
% n = 3 ;
% Si = 2 ;
% Sj = 1 ;

nprey = zeros(3,1) ;
count = 0 ;

% top
z = Si - 1 ;
if z == 0
    z = n ;
end
if M(z,Sj) == 1
    count = count + 1 ;
    nprey(1,count) = 1 ;
    nprey(2,count) = z ;
    nprey(3,count) = Sj ;
end

% right
z = Sj + 1 ;
if z == n + 1
    z = 1 ;
end
if M(Si,z) == 1
    count = count + 1 ;
    nprey(1,count) = 1 ; 
    nprey(2,count) = Si ;
    nprey(3,count) = z ;
end

% bottom
z = Si + 1 ;
if z == n + 1
    z = 1 ;
end
if M(z,Sj) == 1
    count = count + 1 ;
    nprey(1,count) = 1 ;
    nprey(2,count) = z ;
    nprey(3,count) = Sj ;
end
   
% left
z = Sj - 1 ;
if z == 0
    z = n ;
end
if M(Si,z) == 1
    count = count + 1 ;
    nprey(1,count) = 1 ;
    nprey(2,count) = Si ;
    nprey(3,count) = z ;
end

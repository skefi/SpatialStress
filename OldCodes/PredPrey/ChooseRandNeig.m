% Function which chooses randomly one of the 4 neighbors of a prey cell
% Associated with the pg Pascual_AC_WithSkip40.m

function[x] = ChooseRandNeig(Si,Sj)

global n

% Choose a number between 1 and 4
randneig = floor(rand*4) ;
if randneig == 0
    randneig = 4 ;
end

% Associate the number with the coordinate of one of the four neighbors
if randneig == 1 % top
    ni = Si - 1 ;
    if ni == 0
        ni = n ;
    end
    nj = Sj ;

elseif randneig == 2 % right
    ni = Si ;
    nj = Sj + 1 ;
    if nj == n + 1
        nj = 1 ;
    end

elseif randneig == 3 % bottom
    ni = Si + 1 ;
    if ni == n + 1
        ni = 1 ;
    end
    nj = Sj ;

else % Left
    ni = Si ;
    nj = Sj - 1 ;
    if nj == 0
        nj = n ;
    end
end

x = [ni,nj] ; 
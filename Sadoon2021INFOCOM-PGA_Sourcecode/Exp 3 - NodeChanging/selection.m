function [index] = selection(Queen,male,Nm,MQ)

    r = zeros(1,Nm);
    f = zeros(1,Nm);
    F = zeros(1,Nm);
    index = [];

    for j = 1:Nm
        r(j) = pdist2(Queen.Position,male(j).Position);
        f(j) = male(j).Fitness / r(j);
        F(j) = 1 / (1 + f(j));
    end
    P = F(:) ./ sum(F(:));
    Nc = round(rand() * MQ);
    if(Nc >= 1)
        chromosom = cat(1,male.Position);
        index = FuzzyRouletteWheelSelection(chromosom,P,Nc);
    end 

end
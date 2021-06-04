function [index] =  FuzzyRouletteWheelSelection(chromosom,fitness,NumC)

    As = (fitness(:) / sum(fitness)) * 100;
    
    nC = numel(fitness); % Number of Chromosom
    SM = zeros(nC,nC); % Similarity Matrix
    for i = 1:nC
        for j = i:nC
            if(i == j)
                SM(i,j) = nan;
            else
                SM(i,j) = sum(chromosom(i,:) == chromosom(j,:));
                SM(j,i) = SM(i,j);
            end
        end
    end
    
    OW = zeros(5,nC); % Order Matrix: 1 = Chromosom's number, 2 = Value on the plot, 3 = LC, 4 = RC, 5 = CP
    OW(1,1) = 1;
    OW(2,1) = As(1);
    OW(4,1) = As(1);
    OW(5,1) = As(1) / 2;
    G = size(chromosom,2);
    
    FRW = [];
    
    for i = 2:nC
        [IG,OW(1,i)] = max(SM(OW(1,i-1),:));
        SM(OW(1,i-1),OW(1,i)) = nan;
        SM(OW(1,i),OW(1,i-1)) = nan;
        OW(2,i) = OW(2,i-1) + As(OW(1,i));
        OW(3,i) = OW(4,i-1) - (As(OW(1,i-1)) * IG / G);
        OW(4,i) = OW(3,i) + As(OW(1,i));
        OW(5,i) = (OW(3,i) + OW(4,i)) / 2;
        
        FRW(end+1) = (OW(4,i-1) + OW(3,i)) / 2;
    end
    
    r = rand(1,NumC) * OW(4,nC);
    index = zeros(1,NumC);
    if(nC == 1)
        index = OW(1,nC);
    else
        for j = 1:NumC
            if(r(j) > FRW(end))
                index(j) = OW(1,nC);
            else
                for i = 1:size(FRW,2)
                    if(r(j) <= FRW(i))
                        index(j) = OW(1,i);
                        break;                
                    end
                end
            end
        end
    end


end
function [egg] = crossAntMate(queen,male,w,nNode)
    egg = round(abs(queen + rand() * w .* (queen-male)));
    index = find(egg > nNode);
    if(size(index,2) > 0)
        egg(index) = nNode;
    end
    index = find(egg < 1);
    if(size(index,2) > 0)
        egg(index) = 1;
    end
end
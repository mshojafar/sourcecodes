function [egg] = mutateAntMate(egg,M,nNode)
    egg = round(abs(egg + rand() * M .* egg));
    index = find(egg > nNode);
    if(size(index,2) > 0)
        egg(index) = nNode;
    end
%     index = find(egg < 1);
%     if(size(index,2) > 0)
%         egg(index) = 1;
%     end
end
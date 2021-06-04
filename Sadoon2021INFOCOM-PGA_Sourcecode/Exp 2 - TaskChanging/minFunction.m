function [minEng,minComp,minMakespan] = minFunction(task,node)

    minMakespan = sum(task(1,:)) / sum(node(1,:));
    
    minComp = sum(task(1,:)) / max(node(1,:));
    
%     minEng = size(node,2) * min(task(1,:)) * minMakespan;
    minEng = 0;
    for i = 1:size(task,2)
        minEng = min((task(1,i) ./ node(1,:)) .* node(8,:)) + minEng;
    end

end
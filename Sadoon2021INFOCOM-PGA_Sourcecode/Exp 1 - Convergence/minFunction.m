function [minEng,minComp,minMakespan] = minFunction(task,node)

    minMakespan = sum(task(1,:)) / sum(node(1,:));
    
    minComp = sum(task(1,:)) / max(node(1,:));
    
    minEng = size(node,2) * min(node(8,:)) * minMakespan;

end
function [eng] = energyCost(firstAddress,secondAddress,dataSetNum)
    
    node = load('dataset\node.mat');
    activeTime = load(strcat(firstAddress, secondAddress,'.mat'));
    
    makespan = max(activeTime.activeTime,[],2);
    idleTime = repmat(makespan,1,size(node.node , 2)) - activeTime.activeTime;
    
    eng = zeros(dataSetNum,size(node.node , 2));
    
    count = 1;
    while(count <= dataSetNum)
        eng(count,:) = (node.node(8,:) .* activeTime.activeTime(count,:)) + (node.node(7,:) .* idleTime(count,:));
        count = count + 1;
    end
    
    save(strcat(firstAddress,'eng','.mat'),'eng');
    save(strcat(firstAddress,'makespan','.mat'),'makespan');
    

end
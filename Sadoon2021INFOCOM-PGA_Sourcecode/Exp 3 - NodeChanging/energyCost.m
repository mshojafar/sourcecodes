function [eng] = energyCost(firstAddress,secondAddress,dataSetNum)
    
    makespan = zeros(dataSetNum,1);
    counter = 1;
    while(counter <= dataSetNum)
        
        node = load(strcat('dataset\node', num2str(counter),'.mat'));
        activeTime = load(strcat(firstAddress, secondAddress,num2str(counter),'.mat'));
        
        makespan(counter,1) = max(activeTime.activeTime,[],2);
        idleTime = repmat(makespan(counter,1),1,size(node.node , 2)) - activeTime.activeTime;
        
        eng = zeros(1,size(node.node , 2));
        
        eng = (node.node(8,:) .* activeTime.activeTime) + (node.node(7,:) .* idleTime);
        
        save(strcat(firstAddress,'eng',num2str(counter),'.mat'),'eng');
        
        counter = counter + 1;
        
    end 
    
    save(strcat(firstAddress,'makespan','.mat'),'makespan');
    

end
function [] = computationCost(firstAddress,secondAddress,dataSetNum)

    address = strcat(firstAddress, secondAddress);

%     x = load('dataset\dataset.mat');
    
    counter = 1;
    while(counter <= dataSetNum)
        
%         y = load(strcat(address,num2str(counter),'.mat'));
%         node = load(strcat('dataset\node', num2str(counter),'.mat'));
        activeTime = load(strcat(address,'activeTime',num2str(counter),'.mat'));
        
        comp = sum(activeTime.activeTime,2);
        
%         for i = 1:size(y.list , 2)
%             
%             exeTime = x.X(1,i) / node.node(1,y.list(1,i));
%             comp(y.list(1,i)) = comp(y.list(1,i)) + exeTime; 
%         end
        
        save(strcat(address,'comp',num2str(counter),'.mat'),'comp');
        
        counter = counter + 1;
    end 


end
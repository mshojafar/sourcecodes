function MinCalculat(Address,dataSetNum)
    
    x = load('dataset\dataset.mat');
    minEng = zeros(dataSetNum,1);
    minComp = zeros(dataSetNum,1);
    minMakespan = zeros(dataSetNum,1);
    
    counter = 1;
    while(counter <= dataSetNum)
        node = load(strcat('dataset\node', num2str(counter), '.mat'));
        
        [minEng(counter),minComp(counter),minMakespan(counter)] = minFunction(x.X,node.node);
        
        counter = counter + 1;
    end
    
    save(strcat(Address,'minEng','.mat'),'minEng');
    save(strcat(Address,'minComp','.mat'),'minComp');
    save(strcat(Address,'minMakespan','.mat'),'minMakespan');

end
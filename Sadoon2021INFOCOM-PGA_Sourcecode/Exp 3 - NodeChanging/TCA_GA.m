function [gf] = TCA_GA(firstAddress,secondAddress,dataSetNum,MaxRep,nPop,fogNum,cloudNum)
% TCA : Task Classification Algorithm

    x = load('dataset\dataset.mat');
    completeAddress = strcat(firstAddress, secondAddress);
    PDST = zeros(dataSetNum,1);
    
    numT = 1:size(x.X , 2);
    x.X = [numT; x.X(1:end,:)];
    sortData = sortrows(x.X',6)';
    
    DIRatio = zeros(3,size(x.X , 2)); % Deadline Instraction Ratio
    DIRatio(1,:) = sortData(6,:) / max(sortData(6,:)); % Deadline Ratio Row(y)
    DIRatio(2,:) = sortData(2,:) / max(sortData(2,:));% Instraction Ration Row(x)
    
    counter = 1;
    while(counter <= dataSetNum)
        node = load(strcat('dataset\node', num2str(counter), '.mat'));
        list = zeros(2,size(x.X , 2));
        activeTime = zeros(1,size(node.node , 2));
        
        fogList = [];
        cloudList = [];
        
        for i = 1 : size(x.X , 2)
            y = 1 - DIRatio(2,i);
            if(DIRatio(1,i) <= 0.5) && (DIRatio(2,i) <= 0.5)
                DIRatio(3,i) = 1; % 1:fog node , 2:cloud node
                fogList(:,end+1) = sortData(:,i);
            elseif(DIRatio(1,i) > y)
                DIRatio(3,i) = 2; % 1:fog node , 2:cloud node
                cloudList(:,end+1) = sortData(:,i);
            else
                randF = randi([1 fogNum(counter)],1,1);
                randC = randi([fogNum(counter)+1 fogNum(counter)+cloudNum(counter)],1,1);
                e1 = (sortData(2,i) / node.node(1,randF)) * node.node(8,randF);
                e2 = (sortData(2,i) / node.node(1,randC)) * node.node(8,randC);
                if(e1 <= e2)
                    DIRatio(3,i) = 1;
                    fogList(:,end+1) = sortData(:,i);
                else
                    DIRatio(3,i) = 2;
                    cloudList(:,end+1) = sortData(:,i);
                end
            end
        end
        
        TFList = zeros(2,size(fogList , 2));
        TCList = zeros(2,size(cloudList , 2));
        
        [TFList,activeTime(1:fogNum(counter))] = geneticAlgorithm(fogList,node.node(:,1:fogNum(counter)),MaxRep,nPop);
        [TCList,activeTime(fogNum(counter)+1:fogNum(counter)+cloudNum(counter))] = ...
            geneticAlgorithm(cloudList,node.node(:,fogNum(counter)+1:fogNum(counter)+cloudNum(counter)),MaxRep,nPop);
        
        countF = 1;
        countC = 1;
        
        for i = 1:size(x.X , 2)
            if(DIRatio(3,i) == 1)
                list(:,sortData(1,i)) = TFList(:,countF);
                countF = countF + 1;
            else
                list(:,sortData(1,i)) = TCList(:,countC);
                countC = countC + 1;
            end
        end
        
        save(strcat(completeAddress,'activeTime',num2str(counter),'.mat'),'activeTime');
        save(strcat(completeAddress,num2str(counter),'.mat'),'list');
        
        % PDST
        pdstCounter = 0;
        for i = 1:size(x.X , 2)
            if (x.X(6,i) >= list(2,i))
                pdstCounter = pdstCounter + 1;
            end
        end
        PDST(counter,1) = pdstCounter / size(x.X , 2);
        
        counter = counter + 1;
        
    end
    
    save(strcat(completeAddress,'PDST','.mat'),'PDST');
    
    computationCost(firstAddress,secondAddress,dataSetNum);   % Computation and Communication Cost
    energyCost(completeAddress,'activeTime',dataSetNum);
    gf = goalFunction(completeAddress);


end
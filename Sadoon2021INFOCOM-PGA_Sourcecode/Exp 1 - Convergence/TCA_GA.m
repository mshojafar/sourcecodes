function [FogBCH,CloudBCH] = TCA_GA(firstAddress,secondAddress,dataSetNum,MaxRep,nPop,fogNum,cloudNum)
% TCA : Task Classification Algorithm

    node = load('dataset\node.mat');
%     completeAddress = strcat(firstAddress, secondAddress);
    activeTime = zeros(dataSetNum,size(node.node , 2));
    
    counter = 1;
    while(counter <= dataSetNum)
        x = load(strcat('dataset\dataset', num2str(counter), '.mat'));
        list = zeros(2,size(x.X , 2));
        
        numT = [1:size(x.X , 2)];
        x.X = [numT; x.X(1:end,:)];
        sortData = sortrows(x.X',6)';
        
        DIRatio = zeros(3,size(x.X , 2)); % Deadline Instraction Ratio
        DIRatio(1,:) = sortData(6,:) / max(sortData(6,:)); % Deadline Ratio Row(y)
        DIRatio(2,:) = sortData(2,:) / max(sortData(2,:));% Instraction Ration Row(x)
        
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
                randF = randi([1 fogNum],1,1);
                randC = randi([fogNum+1 fogNum+cloudNum],1,1);
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
        
        
        
        [activeTime(counter,1:fogNum),FogBCH] = geneticAlgorithm(fogList,node.node(:,1:fogNum),MaxRep,nPop);
        [activeTime(counter,fogNum+1:fogNum+cloudNum),CloudBCH] = ...
            geneticAlgorithm(cloudList,node.node(:,fogNum+1:fogNum+cloudNum),MaxRep,nPop);
        
        counter = counter + 1;
        
    end


end
function [gf] = power2choices(firstAddress,secondAddress,dataSetNum)
    
    x = load('dataset\dataset.mat');
    completeAddress = strcat(firstAddress, secondAddress);
    PDST = zeros(dataSetNum,1);
    
    counter = 1;
    while(counter <= dataSetNum)
        
        node = load(strcat('dataset\node', num2str(counter), '.mat'));
        list = zeros(2,size(x.X , 2));
        activeTime = zeros(1,size(node.node , 2));
        
        ctrl = zeros(3,2);
        ctrlTotal = zeros(3,size(node.node , 2));
        
        for i = 1:size(x.X , 2)
            
            sizeA = size(node.node , 2);
            twoChoice = randi([1 sizeA],1,2);
            
            while(twoChoice(1) == twoChoice(2))
                twoChoice(1) = randi(sizeA);
            end
            ctrl(:,1) = ctrlTotal(:,twoChoice(1));
            ctrl(:,2) = ctrlTotal(:,twoChoice(2));
            
            for j = 1:2
                    
                exeTime = x.X(1,i) / node.node(1,twoChoice(j));
%                 compNode = (exeTime * node.node(2,twoChoice(j))) + (x.X(2,i) * node.node(3,twoChoice(j)));
%                 compCost = compNode;
                
                resTime = (2 * node.node(9,twoChoice(j))) + exeTime + activeTime(twoChoice(j));
                
%                 top = resTime - x.X(5,i);
%                 if(top > 0)
%                     violation = (top / x.X(5,i)) * 100;
%                 end
%                 violCost = x.X(7,i) * violation;
                
                ctrl(1,j) = twoChoice(j);
%                 ctrl(2,j) = compCost;
                ctrl(2,j) = ctrl(2,j) + exeTime;
%                 ctrl(4,j) = violCost;
                ctrl(3,j) = resTime;
                
%                 violCost = 0;
%                 violation = 0;
                
            end
            
            if(ctrl(3,1) <= ctrl(3,2))
                list(1,i) = ctrl(1,1);
                list(2,i) = ctrl(3,1);
                activeTime(ctrl(1,1)) = ctrl(2,1);
                ctrlTotal(1:3,ctrl(1,1)) = ctrl(1:3,1);
            else
                list(1,i) = ctrl(1,2);
                list(2,i) = ctrl(3,2);
                activeTime(ctrl(1,2)) = ctrl(2,2);
                ctrlTotal(1:3,ctrl(1,2)) = ctrl(1:3,2);
            end
            
            ctrl(:) = 0;
            
        end
        
        save(strcat(completeAddress,'activeTime',num2str(counter),'.mat'),'activeTime');
        save(strcat(completeAddress,num2str(counter),'.mat'),'list');
        
        
        ctrlTotal(:) = 0;
        
        % PDST
        pdstCounter = 0;
        for i = 1:size(x.X , 2)
            if (x.X(5,i) >= list(2,i))
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
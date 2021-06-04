function [gf] = power2choices(firstAddress,secondAddress,dataSetNum)

    node = load('dataset\node.mat');
    completeAddress = strcat(firstAddress, secondAddress);
    activeTime = zeros(dataSetNum,size(node.node , 2));
    PDST = zeros(dataSetNum,1);
    
    ctrl = zeros(3,2);
    ctrlTotal = zeros(3,size(node.node , 2));
    
    counter = 1;
    while(counter <= dataSetNum)
        
        x = load(strcat('dataset\dataset', num2str(counter), '.mat'));
        list = zeros(2,size(x.X , 2));
        
        for i = 1:size(x.X , 2)
            
            sizeA = size(node.node , 2);
            twoChoice = randi([1 sizeA],1,2);
            
            ctrl(:,1) = ctrlTotal(:,twoChoice(1));
            ctrl(:,2) = ctrlTotal(:,twoChoice(2));
            
            for j = 1:2
                    
                exeTime = x.X(1,i) / node.node(1,twoChoice(j));
                
                resTime = (2 * node.node(9,twoChoice(j))) + exeTime + activeTime(counter,twoChoice(j));
                
                ctrl(1,j) = twoChoice(j);
                ctrl(2,j) = ctrl(2,j) + exeTime;
                ctrl(3,j) = resTime;
                
            end
            
            if(ctrl(3,1) <= ctrl(3,2))
                list(1,i) = ctrl(1,1);
                list(2,i) = ctrl(3,1);
                activeTime(counter,ctrl(1,1)) = ctrl(2,1);
                ctrlTotal(1:3,ctrl(1,1)) = ctrl(1:3,1);
            else
                list(1,i) = ctrl(1,2);
                list(2,i) = ctrl(3,2);
                activeTime(counter,ctrl(1,2)) = ctrl(2,2);
                ctrlTotal(1:3,ctrl(1,2)) = ctrl(1:3,2);
            end
            
            ctrl(:) = 0;
            
        end
        
        data_name = strcat(completeAddress,num2str(counter),'.mat');
        save(data_name,'list');
        
        
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
    
    save(strcat(completeAddress,'activeTime','.mat'),'activeTime');
    save(strcat(completeAddress,'PDST','.mat'),'PDST');
    
    computationCost(firstAddress,secondAddress);   % Computation and Communication Cost
    energyCost(completeAddress,'activeTime',dataSetNum);
    gf = goalFunction(completeAddress);
            

end
function [activeTime,BCH] = geneticAlgorithm(task,node,MaxRep,nPop)

       
    nNode = size(node , 2);
    nTask = size(task , 2);
    
    activeTime = zeros(1,nNode);  
    list = zeros(2,nTask);
    
    empty_Gen.Position = [];
    empty_Gen.Fitness = inf;
    empty_Gen.Comp = inf;
    empty_Gen.Eng = inf;
    empty_Gen.ActiveTime = [];
    empty_Gen.Response = [];

    Gen = repmat(empty_Gen,nPop,1);
    
    ChangedTask = task;
    ChangedTask(1,:) = [];
    [minEng,minComp,~] = minFunction(ChangedTask,node);
    
    for i = 1:nPop
        % Initialize Position
        Gen(i).Position = randi([1,nNode],1,nTask);

        % Evaluate
        [Gen(i).Fitness,Gen(i).Comp,Gen(i).Eng,Gen(i).ActiveTime,Gen(i).Response] = ...
            costFunction(Gen(i).Position,task,node,minEng,minComp);
    end
    
    [bestAns,IPop] = max([Gen.Fitness]);
    list(1,:) = Gen(IPop).Position;
    list(2,:) = Gen(IPop).Response;
    activeTime = Gen(IPop).ActiveTime;
    
    BCH = zeros(1,100);
    indexB = 1;
    
    %Number of Generations
    for i = 1:MaxRep

        Gen = maskCrossOver(Gen,task,node,minEng,minComp,nPop);
        Gen = mutation(Gen,task,node,minEng,minComp);
        
        if(bestAns < max([Gen.Fitness]))
            [bestAns,IPop] = max([Gen.Fitness]);
            list(1,:) = Gen(IPop).Position;
            list(2,:) = Gen(IPop).Response;
            activeTime = Gen(IPop).ActiveTime;
        end
        
        if((i/10) == indexB)
            BCH(indexB) = bestAns;
            indexB = indexB + 1;
        end

    end
    
%     figure
%     plot(1:MaxRep,BCH);
%     xlabel('Iteration');
%     ylabel('Fitness');
%     title(strcat('Convergence:',layer,' , ',num2str(nPop)));
%     
%     save(strcat(address,layer),'BCH');
    
end
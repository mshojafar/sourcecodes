function [gf] = AMO(firstAddress,secondAddress,dataSetNum,MaxRep,nPop)

    node = load('dataset\node.mat');
    nNode = size(node.node , 2);
    completeAddress = strcat(firstAddress, secondAddress);
    activeTime = zeros(dataSetNum,nNode);
    PDST = zeros(dataSetNum,1);
    
    minEng = load('MIterations\minEng.mat');
    minMakespan = load('MIterations\minMakespan.mat');
    
    empty_Ant.Position = [];
    empty_Ant.Fitness = inf;
    empty_Ant.Makespan = inf;
    empty_Ant.Eng = inf;
    empty_Ant.ActiveTime = [];
    empty_Ant.Response = [];

    counter = 1;
    while(counter <= dataSetNum)
        x = load(strcat('dataset\dataset', num2str(counter), '.mat'));
        nTask = size(x.X , 2);
        list = zeros(2,nTask);
        
        Ant = repmat(empty_Ant,nPop,1);

        for i = 1:nPop
            % Initialize Position
            Ant(i).Position = randi([1,nNode],1,nTask);

            % Evaluate
            [Ant(i).Fitness,Ant(i).Eng,Ant(i).Makespan,Ant(i).ActiveTime,Ant(i).Response] = ...
                AntFitness(Ant(i).Position,x.X,node.node,minEng.minEng(counter),minMakespan.minMakespan(counter));
        end
        
        PQ = (0.4-0.2).*rand() + 0.2; % Percent of Queen 
        w = rand();
        EF = randi([nPop nPop+200]); % Extinction Factor 
        M = 0.05;
        OF = 0.9;
        cmf = 0.2;

        for k = 1:MaxRep 
            Nq = round(PQ * numel(Ant));
            Nm = numel(Ant) - Nq;
            MQ = round(Nm * rand());
%             Queen = repmat(empty_Ant,Nq,1);
%             male = repmat(empty_Ant,Nm,1);
            egg = repmat(empty_Ant,1);
            [~,I] = maxk([Ant.Fitness],Nq);
            Queen = Ant(I);
            IN = 1:numel(Ant);
            IN(I(:)) = [];
            male = Ant(IN);
            
            for i = 1:Nq
                index = selection(Queen(i),male,Nm,MQ);
                index(index == 0) = [];
                if(numel(index) ~= 0)
                    for j = 1:numel(index)
                        r = rand();
                        if(OF > r)
                            egg(end+1).Position = crossAntMate([Queen(i).Position],[male(index(j)).Position],w,nNode);
                            if(cmf > r)
                                egg(end).Position = mutateAntMate([egg(end).Position],M,nNode);    
                            end
                            [egg(end).Fitness,egg(end).Eng,egg(end).Makespan,egg(end).ActiveTime,egg(end).Response] = ...
                                    AntFitness(egg(end).Position,x.X,node.node,minEng.minEng(counter),minMakespan.minMakespan(counter));
                        end
                    end 
                end
            end
            egg(1) = [];
            Ant(:) = [];
            Ant(1:numel(Queen)) = Queen; 
            col = numel(Queen)+numel(male);
            Ant(end+1:col) = male;
            Ant(end+1:col+numel(egg)) = egg;
            if(numel(Ant) > EF)
                [~,I] = maxk([Ant.Fitness],EF);
                Ant = Ant(I);
            end
            Queen(:) = [];
            male(:) = [];
        end
        
        [~,I] = max([Ant.Fitness]);
        list(1,:) = Ant(I).Position;
        list(2,:) = Ant(I).Response;
        
        data_name = strcat(completeAddress,num2str(counter),'.mat');
        save(data_name,'list');
        
        activeTime(counter,:) = [Ant(I).ActiveTime];
        
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
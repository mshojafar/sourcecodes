function [Gen] = maskCrossOver(Gen,task,node,minEng,minComp,nPop)

    percent = (1-0.8) * rand() + 0.8;
    r = round(percent * nPop / 2);
    randID = randi([1 nPop],1,r);
    mask = randi([0 1],r,size(task,2));
    child1 = zeros(1,size(task,2));
    child2 = zeros(1,size(task,2));
    
    for i = 1:r
        index = RouletteWheelSelection([Gen.Fitness]);
        
        for j = 1:size(task,2)
            if(mask(i,j) == 0)
                child1(j) = Gen(index).Position(j);
                child2(j) = Gen(randID(i)).Position(j);
            elseif(mask(i,j) == 1)
                child1(j) = Gen(randID(i)).Position(j);
                child2(j) = Gen(index).Position(j);
            end
        end
        
        Gen(end+1).Position = child1;
        [Gen(end).Fitness,Gen(end).Comp,Gen(end).Eng,Gen(end).ActiveTime,Gen(end).Response] = costFunction(Gen(end).Position,task,node,minEng,minComp);
        Gen(end+1).Position = child2;
        [Gen(end).Fitness,Gen(end).Comp,Gen(end).Eng,Gen(end).ActiveTime,Gen(end).Response] = costFunction(Gen(end).Position,task,node,minEng,minComp);
        
    end
    
    [~,IPop] = maxk([Gen.Fitness],nPop);

    Gen = Gen(IPop(:));

end
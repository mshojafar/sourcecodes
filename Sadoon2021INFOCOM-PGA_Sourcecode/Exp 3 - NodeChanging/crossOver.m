function [Gen] = crossOver(Gen,task,node,minEng,minComp,nPop)

    r = randi([150 nPop]);
    randID = randi([1 nPop],1,r);
    
    for i = 1:r
        
        index = RouletteWheelSelection([Gen.Fitness]);
        
        i1 = randi([1 (size(task , 2) - 2)]);
        i2 = randi([i1 (size(task , 2) - 1)]);
        
        Gen(end+1).Position = [Gen(index).Position(1:i1) Gen(randID(i)).Position(i1+1:i2) Gen(index).Position(i2+1:size(task , 2))];
        [Gen(end).Fitness,Gen(end).Comp,Gen(end).Eng,Gen(end).ActiveTime,Gen(end).Response] = costFunction(Gen(end).Position,task,node,minEng,minComp);
        
        Gen(end+1).Position = [Gen(randID(i)).Position(1:i1) Gen(index).Position(i1+1:i2) Gen(randID(i)).Position(i2+1:size(task , 2))];
        [Gen(end).Fitness,Gen(end).Comp,Gen(end).Eng,Gen(end).ActiveTime,Gen(end).Response] = costFunction(Gen(end).Position,task,node,minEng,minComp);
    end
    
    [~,IPop] = maxk([Gen.Fitness],nPop);

    Gen = Gen(IPop(:));


end
function [Gen] = mutation(Gen,task,node,minEng,minComp)
    
    numberKromozom = round(size(Gen , 1) * 0.02);
    for i = 1:numberKromozom

        r = randi([1 size(Gen , 1)]);
        randPosition = randi([1 size(task , 2)]);
        randNode = randi([1 size(node , 2)]);
        Gen(r).Position(randPosition) = randNode;
        [Gen(r).Fitness,Gen(r).Comp,Gen(r).Eng,Gen(r).ActiveTime,Gen(r).Response] = costFunction(Gen(r).Position,task,node,minEng,minComp);
    end

end
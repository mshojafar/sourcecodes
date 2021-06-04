function [Fitness,Comp,Eng,ActiveTime,Response] = costFunction(position,task,node,minEng,minComp)
    
    ActiveTime = zeros(1,size(node,2));
    Response = zeros(1,size(task,2));
    Comp = 0;
    pdstCounter = 0;
    
    for i = 1:size(task , 2)
            
        exeTime = task(2,i) / node(1,position(i));
        Comp = Comp + exeTime;
        ActiveTime(position(i)) = ActiveTime(position(i)) + exeTime;
        Response(i) = ActiveTime(position(i)) + (2 * node(9,position(i)));
        
        if (task(6,i) >= Response(i))
            pdstCounter = pdstCounter + 1;
        end
        
    end
    
    PDST = pdstCounter / size(task , 2);
    
    makespan = max(ActiveTime);
    IdleTime = makespan - ActiveTime(:);
    Eng = sum((node(8,:) .* ActiveTime(:)) + (node(7,:) .* IdleTime(:)));
    
    
    Fitness = 0.33 * (minEng / sum(Eng)) + 0.33 * (minComp / Comp) + 0.33 * PDST;
    


end
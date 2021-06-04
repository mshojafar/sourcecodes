function [Fitness,Eng,Makespan,ActiveTime,Response] = AntFitness(position,task,node,minEng,minMakespan)

    
    ActiveTime = zeros(1,size(node,2));
    Response = zeros(1,size(task,2));
    Comp = 0;
    
    for i = 1:size(task , 2)
            
        exeTime = task(1,i) / node(1,position(i));
        Comp = Comp + exeTime;
        ActiveTime(position(i)) = ActiveTime(position(i)) + exeTime;
        Response(i) = ActiveTime(position(i)) + (2 * node(9,position(i)));
        
    end
    
    Makespan = max(ActiveTime);
    IdleTime = Makespan - ActiveTime(:);
    Eng = sum((node(8,:) .* ActiveTime(:)) + (node(7,:) .* IdleTime(:)));
    
    
    Fitness = 0.5 * (minEng / sum(Eng)) + (1 - 0.5) * (minMakespan / Makespan);

end
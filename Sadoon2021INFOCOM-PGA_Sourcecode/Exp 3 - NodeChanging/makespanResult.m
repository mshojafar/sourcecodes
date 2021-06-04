function [] = makespanResult(fAP2C,fAANT,fATCAGA,iteration,dataSetNum)

    rMP2C = zeros(iteration,dataSetNum);
    rMANT = zeros(iteration,dataSetNum);
    rMTCAGA = zeros(iteration,dataSetNum);

    counter = 1;
    while(counter <= iteration)
        
        secondAddress = strcat(num2str(counter),'\','makespan.mat');
        
        MP2C = load(strcat(fAP2C,secondAddress));
        MANT = load(strcat(fAANT,secondAddress));
        MTCAGA = load(strcat(fATCAGA,secondAddress));

        rMP2C(counter,:) = MP2C.makespan;
        rMANT(counter,:) = MANT.makespan;
        rMTCAGA(counter,:) = MTCAGA.makespan;
        
        counter = counter + 1;
    end
    
    ctrl = zeros(3,dataSetNum);
    ctrl(1,:) = sum(rMP2C,1) / iteration;
    ctrl(2,:) = sum(rMANT,1) / iteration;
    ctrl(3,:) = sum(rMTCAGA,1) / iteration;
    
    save(strcat('MIterations\','MakespanResult.mat'),'ctrl');
    
    x = [20 40 60 80 100];

    vals = [ctrl(1,1) ctrl(2,1) ctrl(3,1);
            ctrl(1,2) ctrl(2,2) ctrl(3,2);
            ctrl(1,3) ctrl(2,3) ctrl(3,3);
            ctrl(1,4) ctrl(2,4) ctrl(3,4);
            ctrl(1,5) ctrl(2,5) ctrl(3,5);];
        
        
    figure
    bar(x,vals);
    legend('Power of 2 Choices','ANT','TCAGA');
    title('Makespan Result');


end
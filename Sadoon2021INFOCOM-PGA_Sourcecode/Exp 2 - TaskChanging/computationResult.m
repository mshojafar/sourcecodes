function [] = computationResult(fAP2C,fAANT,fATCAGA,iteration,dataSetNum)
    
    rCompP2C = zeros(iteration,dataSetNum);
    rCompANT = zeros(iteration,dataSetNum);
    rCompTCAGA = zeros(iteration,dataSetNum);

    counter = 1;
    while(counter <= iteration)
        
        secondAddress = strcat(num2str(counter),'\','comp.mat');
        
        compP2C = load(strcat(fAP2C,secondAddress));
        compANT = load(strcat(fAANT,secondAddress));
        compTCAGA = load(strcat(fATCAGA,secondAddress));

        rCompP2C(counter,:) = sum(compP2C.comp,2);
        rCompANT(counter,:) = sum(compANT.comp,2);
        rCompTCAGA(counter,:) = sum(compTCAGA.comp,2);
        
        counter = counter + 1;
    end
    
    ctrl = zeros(3,dataSetNum);
    ctrl(1,:) = sum(rCompP2C,1) / iteration;
    ctrl(2,:) = sum(rCompANT,1) / iteration;
    ctrl(3,:) = sum(rCompTCAGA,1) / iteration;
    
    save(strcat('MIterations\','ComputationResult.mat'),'ctrl');
    
    x = [100 200 300 400 500];

    vals = [ctrl(1,1) ctrl(2,1) ctrl(3,1);
            ctrl(1,2) ctrl(2,2) ctrl(3,2);
            ctrl(1,3) ctrl(2,3) ctrl(3,3);
            ctrl(1,4) ctrl(2,4) ctrl(3,4);
            ctrl(1,5) ctrl(2,5) ctrl(3,5);];
        
        
    figure
    bar(x,vals);
    legend('Power of 2 Choices','ANT-mate','TCAGA');
    title('Computation Result');

end
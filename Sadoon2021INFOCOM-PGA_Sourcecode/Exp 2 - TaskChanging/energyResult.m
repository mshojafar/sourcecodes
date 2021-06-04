function [] = energyResult(fAP2C,fAANT,fATCAGA,iteration,dataSetNum)
    
    rEngP2C = zeros(iteration,dataSetNum);
    rEngANT = zeros(iteration,dataSetNum);
    rEngTCAGA = zeros(iteration,dataSetNum);

    counter = 1;
    while(counter <= iteration)
        
        secondAddress = strcat(num2str(counter),'\','eng.mat');
        
        engP2C = load(strcat(fAP2C,secondAddress));
        engANT = load(strcat(fAANT,secondAddress));
        engTCAGA = load(strcat(fATCAGA,secondAddress));
  
        rEngP2C(counter,:) = sum(engP2C.eng,2);
        rEngANT(counter,:) = sum(engANT.eng,2);
        rEngTCAGA(counter,:) = sum(engTCAGA.eng,2);
        
        counter = counter + 1;
    end
    
    ctrl = zeros(3,dataSetNum);
    ctrl(1,:) = sum(rEngP2C,1) / iteration;
    ctrl(2,:) = sum(rEngANT,1) / iteration;
    ctrl(3,:) = sum(rEngTCAGA,1) / iteration;
    
    save(strcat('MIterations\','EnergyResult.mat'),'ctrl');
    
    x = [100 200 300 400 500];

    vals = [ctrl(1,1) ctrl(2,1) ctrl(3,1);
            ctrl(1,2) ctrl(2,2) ctrl(3,2);
            ctrl(1,3) ctrl(2,3) ctrl(3,3);
            ctrl(1,4) ctrl(2,4) ctrl(3,4);
            ctrl(1,5) ctrl(2,5) ctrl(3,5);];
           
    figure
    bar(x,vals);
    legend('Power of 2 Choices','ANT-mate','TCAGA');
    title('Energy Result');


end
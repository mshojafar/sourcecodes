function [] = main()

    fAP2C = 'MIterations\power2choices\';
    fAANT = 'MIterations\antmate\';
    fATCAGA = 'MIterations\tcaga\';

    iteration = 5;
    MaxRep = 50; % Max repitation of meta algorithms
    nPop = 400;
    dataTask = 250; % [20 40 80 160 320]
    fogNum = [15 30 45 60 75];
    cloudNum = [5 10 15 20 25];
    dataSetNum = numel(fogNum);
    
    resultP2C = zeros(iteration,dataSetNum);
    resultANT = zeros(iteration,dataSetNum);
    resultTCAGA = zeros(iteration,dataSetNum);
    

    counter = 1;
    while(counter <= iteration)
        
        secondAddress = strcat(num2str(counter),'\');
        
        adrsP2C = strcat(fAP2C,secondAddress);
        adrsANT = strcat(fAANT,secondAddress);
        adrsTCAGA = strcat(fATCAGA,secondAddress);
        
        status = mkdir(adrsP2C);
        status = mkdir(adrsANT);
        status = mkdir(adrsTCAGA);
   
        createDataset(dataTask);
        createNode(fogNum,cloudNum,dataSetNum);
        MinCalculat('MIterations\',dataSetNum);

        gfP2C = power2choices(fAP2C,secondAddress,dataSetNum);
        gfANT = AMO(fAANT,secondAddress,dataSetNum,MaxRep,nPop);
        gfTCAGA = TCA_GA(fATCAGA,secondAddress,dataSetNum,MaxRep,nPop,fogNum,cloudNum);
        
        resultP2C(counter,:) = gfP2C;
        resultANT(counter,:) = gfANT;
        resultTCAGA(counter,:) = gfTCAGA;
        
        counter = counter + 1;
    end
    
    save(strcat(fAP2C,'GFResult','.mat'),'resultP2C');
    save(strcat(fAANT,'GFResult','.mat'),'resultANT');
    save(strcat(fATCAGA,'GFResult','.mat'),'resultTCAGA');
        
    ctrl = zeros(3,dataSetNum);
    ctrl(1,:) = sum(resultP2C,1) / iteration;
    ctrl(2,:) = sum(resultANT,1) / iteration;
    ctrl(3,:) = sum(resultTCAGA,1) / iteration;
    
    save(strcat('MIterations\','GFResult.mat'),'ctrl');
    
    x = [20 40 60 80 100];

    vals = [ctrl(1,1) ctrl(2,1) ctrl(3,1);
            ctrl(1,2) ctrl(2,2) ctrl(3,2);
            ctrl(1,3) ctrl(2,3) ctrl(3,3);
            ctrl(1,4) ctrl(2,4) ctrl(3,4);
            ctrl(1,5) ctrl(2,5) ctrl(3,5);];
        
        
    figure
    bar(x,vals);
    legend('Power of 2 Choices','ANT-mate','TCAGA');
    title('Goal Function Result');

    computationResult(fAP2C,fAANT,fATCAGA,iteration,dataSetNum);
    makespanResult(fAP2C,fAANT,fATCAGA,iteration,dataSetNum);
    PDST(fAP2C,fAANT,fATCAGA,iteration,dataSetNum);
    energyResult(fAP2C,fAANT,fATCAGA,iteration,dataSetNum);


end
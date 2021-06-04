function [] = main()

    fATCAGA = 'Plots\';

    RUN = 1;
    MaxRep = 10; % Max repitation of meta algorithms
    nPop = [50 100 200 400];
    dataTask = 250; % [20 40 80 160 320]
    fogNum = 45;
    cloudNum = 15;
    dataSetNum = numel(dataTask);
    
    FogBCH = zeros(RUN,100);
    CloudBCH = zeros(RUN,100);
    
    counter = 1;
    while(counter <= numel(nPop))
        
        secondAddress = strcat(num2str(counter),'\');
        adrsTCAGA = strcat(fATCAGA,secondAddress);
        status = mkdir(adrsTCAGA);
        
        for i = 1:RUN
            createDataset(dataTask,dataSetNum);
            createNode(fogNum,cloudNum);
            [FogBCH(i,:),CloudBCH(i,:)] = TCA_GA(fATCAGA,secondAddress,dataSetNum,MaxRep,nPop(counter),fogNum,cloudNum); 
        end
        
        FogBCH = sum(FogBCH,1)/ RUN;
        CloudBCH = sum(CloudBCH,1)/ RUN;
        
        save(strcat(adrsTCAGA,'FogBCH'),'FogBCH');
        save(strcat(adrsTCAGA,'CloudBCH'),'CloudBCH');
        
        counter = counter + 1;
    end


end
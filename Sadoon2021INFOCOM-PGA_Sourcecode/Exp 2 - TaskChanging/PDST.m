function [] = PDST(fAP2C,fAANT,fATCAGA,iteration,dataSetNum)

    
    rPDSTP2C = zeros(iteration,dataSetNum);
    rPDSTANT = zeros(iteration,dataSetNum);
    rPDSTTCAGA = zeros(iteration,dataSetNum);
    
    counter = 1;
    while(counter <= iteration)
        
        secondAddress = strcat(num2str(counter),'\','PDST.mat');
        
        PDSTP2C = load(strcat(fAP2C,secondAddress));
        PDSTANT = load(strcat(fAANT,secondAddress));
        PDSTTCAGA = load(strcat(fATCAGA,secondAddress));
        
        rPDSTP2C(counter,:) = PDSTP2C.PDST;
        rPDSTANT(counter,:) = PDSTANT.PDST;
        rPDSTTCAGA(counter,:) = PDSTTCAGA.PDST;
        
        counter = counter + 1;
        
    end
    
    ctrl = zeros(3,dataSetNum);
    ctrl(1,:) = sum(rPDSTP2C,1) / iteration;
    ctrl(2,:) = sum(rPDSTANT,1) / iteration;
    ctrl(3,:) = sum(rPDSTTCAGA,1) / iteration;
    
    save(strcat('MIterations\','PDSTResult.mat'),'ctrl'); % Each row related to specific method.
    
    x = [100 200 300 400 500];
    
    vals = [ctrl(1,1) ctrl(2,1) ctrl(3,1);
            ctrl(1,2) ctrl(2,2) ctrl(3,2);
            ctrl(1,3) ctrl(2,3) ctrl(3,3);
            ctrl(1,4) ctrl(2,4) ctrl(3,4);
            ctrl(1,5) ctrl(2,5) ctrl(3,5);];
        
        
    figure
    bar(x,vals);
    legend('Power of 2 Choices','ANT-mate','TCAGA');
    title('PDST');


end
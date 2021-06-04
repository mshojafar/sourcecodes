function createNode(fogNum,cloudNum)

    N = fogNum + cloudNum;
    node = zeros(10,N);
    

    %CPU clock rate
    node(1,1:fogNum) = randi([500 2000],1,fogNum);
    node(1,fogNum+1:N) = randi([3000 10000],1,cloudNum);

    %CPU usage cost
    node(2,1:fogNum) = (0.4-0.1).*rand(1,fogNum) + 0.1;
    node(2,fogNum+1:N) = (1-0.7).*rand(1,cloudNum) + 0.7;

    %Memory usage cost
    node(3,1:fogNum) = (0.03-0.01).*rand(1,fogNum) + 0.01;
    node(3,fogNum+1:N) = (0.05-0.02).*rand(1,cloudNum) + 0.02;

    %Bandwidth usage cost
    node(4,1:fogNum) = (0.02-0.01).*rand(1,fogNum) + 0.01;
    node(4,fogNum+1:N) = (0.1-0.05).*rand(1,cloudNum) + 0.05;
    
    %Energy usage cost
    node(5,1:fogNum) = (0.12-0.006).*rand(1,fogNum) + 0.006;
    node(5,fogNum+1:N) = (0.36-0.018).*rand(1,cloudNum) + 0.018;
    
    %Memory
    node(6,1:fogNum) = (250-150).*rand(1,fogNum) + 150;
    node(6,fogNum+1:N) = (2000-200).*rand(1,cloudNum) + 200;
    
    
    %Max Power
    node(8,1:fogNum) = (130-80).*rand(1,fogNum) + 80;
    node(8,fogNum+1:N) = (500-150).*rand(1,cloudNum) + 150;
    
    %Min Power
    node(7,1:fogNum) = (((70-60).*rand(1,fogNum) + 60)/ 100) .* node(8,1:fogNum);
    node(7,fogNum+1:N) = (((70-60).*rand(1,cloudNum) + 60)/ 100) .* node(8,fogNum+1:N);
    
    %delay
    node(9,1:fogNum) = (0.01-0.001).*rand(1,fogNum) + 0.001;
    node(9,fogNum+1:N) = (0.5-0.1).*rand + 0.1;
    
    %Sleep Power
    node(10,1:fogNum) = (((5-1).*rand(1,fogNum) + 1)/ 100) .* node(8,1:fogNum);
    node(10,fogNum+1:N) = (((5-1).*rand(1,cloudNum) + 1)/ 100) .* node(8,fogNum+1:N);


    save('dataset\node.mat','node');


end
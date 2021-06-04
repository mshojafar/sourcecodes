function createNode(fogNum,cloudNum,n)

    N = zeros(1,n);
    
    for i = 1:n
        
        N(i) = fogNum(i) + cloudNum(i);
        node = zeros(10,N(i));
        
        %CPU clock rate
        node(1,1:fogNum(i)) = randi([500 2000],1,fogNum(i));
        node(1,fogNum(i)+1:N(i)) = randi([3000 10000],1,cloudNum(i));

        %CPU usage cost
        node(2,1:fogNum(i)) = (0.4-0.1).*rand(1,fogNum(i)) + 0.1;
        node(2,fogNum(i)+1:N(i)) = (1-0.7).*rand(1,cloudNum(i)) + 0.7;

        %Memory usage cost
        node(3,1:fogNum(i)) = (0.03-0.01).*rand(1,fogNum(i)) + 0.01;
        node(3,fogNum(i)+1:N(i)) = (0.05-0.02).*rand(1,cloudNum(i)) + 0.02;

        %Bandwidth usage cost
        node(4,1:fogNum(i)) = (0.02-0.01).*rand(1,fogNum(i)) + 0.01;
        node(4,fogNum(i)+1:N(i)) = (0.1-0.05).*rand(1,cloudNum(i)) + 0.05;

        %Energy usage cost
        node(5,1:fogNum(i)) = (0.12-0.006).*rand(1,fogNum(i)) + 0.006;
        node(5,fogNum(i)+1:N(i)) = (0.36-0.018).*rand(1,cloudNum(i)) + 0.018;

        %Memory
        node(6,1:fogNum(i)) = (250-150).*rand(1,fogNum(i)) + 150;
        node(6,fogNum(i)+1:N(i)) = (2000-200).*rand(1,cloudNum(i)) + 200;


        %Max Power
        node(8,1:fogNum(i)) = (130-80).*rand(1,fogNum(i)) + 80;
        node(8,fogNum(i)+1:N(i)) = (500-150).*rand(1,cloudNum(i)) + 150;

        %Min Power
        node(7,1:fogNum(i)) = (((70-60).*rand(1,fogNum(i)) + 60)/ 100) .* node(8,1:fogNum(i));
        node(7,fogNum(i)+1:N(i)) = (((70-60).*rand(1,cloudNum(i)) + 60)/ 100) .* node(8,fogNum(i)+1:N(i));

        %delay
        node(9,1:fogNum(i)) = (0.01-0.001).*rand(1,fogNum(i)) + 0.001;
        node(9,fogNum(i)+1:N(i)) = (0.5-0.1).*rand + 0.1;

        %Sleep Power
        node(10,1:fogNum(i)) = (((5-1).*rand(1,fogNum(i)) + 1)/ 100) .* node(8,1:fogNum(i));
        node(10,fogNum(i)+1:N(i)) = (((5-1).*rand(1,cloudNum(i)) + 1)/ 100) .* node(8,fogNum(i)+1:N(i));
        
        save(strcat('dataset\node',num2str(i),'.mat'),'node');
        
    end


end
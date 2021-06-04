function createDataset(n)
    
    
    X = zeros(7,n);
    X(2,:) = (200-50).*rand(1,n) + 50; %Memory required
    X(3,:) = (100-10).*rand(1,n) + 10; %Input file size
    X(4,:) = (100-10).*rand(1,n) + 10; %Output file size
    X(6,:) = (99.99-90).*rand(1,n) + 90; %Quality of Service
    X(7,:) = (5-2).*rand(1,n) + 2; %Penalty

    for j = 1:n
        randNum = randi(3);
        if (randNum == 1)
            X(1,j) = randi([100 372]); %Number of instructions(MI)
            X(5,j) = (0.5-0.1).*rand(1) + 0.1; %Deadline
        elseif (randNum == 2)
            X(1,j) = randi([1028 4280]); %Number of instructions(MI)
            X(5,j) = (2.5-0.5).*rand(1) + 0.5; %Deadline
        elseif (randNum == 3)
            X(1,j) = randi([5123 9784]); %Number of instructions(MI)
            X(5,j) = (10-2.5).*rand(1) + 2.5; %Deadline
        end
    end

    save('dataset\dataset.mat','X');

end
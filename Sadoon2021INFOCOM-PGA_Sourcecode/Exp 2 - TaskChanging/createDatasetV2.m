function createDatasetV2(n,dataSetNum)


    for i = 1:dataSetNum
        X = zeros(7,n(i));
        X(1,:) = randi([100 10000],1,n(i)); %Number of instructions(MI)
        X(2,:) = (200-50).*rand(1,n(i)) + 50; %Memory required
        X(3,:) = (100-10).*rand(1,n(i)) + 10; %Input file size
        X(4,:) = (100-10).*rand(1,n(i)) + 10; %Output file size
        X(5,:) = (10-0.1).*rand(1,n(i)) + 0.1; %Deadline
        X(6,:) = (99.99-90).*rand(1,n(i)) + 90; %Quality of Service
        X(7,:) = (5-2).*rand(1,n(i)) + 2; %Penalty

        data_name = strcat('dataset\dataset',int2str(i),'.mat');
        save(data_name,'X');
    end

end
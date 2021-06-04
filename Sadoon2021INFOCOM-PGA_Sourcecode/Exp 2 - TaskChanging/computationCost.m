function [] = computationCost(firstAddress,secondAddress)

    address = strcat(firstAddress, secondAddress);    
    activeTime = load(strcat(address,'activeTime.mat'));
    comp = sum(activeTime.activeTime,2);
    save(strcat(address,'comp','.mat'),'comp');
end
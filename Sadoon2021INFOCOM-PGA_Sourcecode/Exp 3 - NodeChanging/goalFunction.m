function [F] = goalFunction(firstAddress)

    PDST = load(strcat(firstAddress,'PDST','.mat'));
    minEng = load('MIterations\minEng.mat');
    minComp = load('MIterations\minComp.mat');
    
    F = zeros(1,5);
    
    for i = 1:5
        computation = load(strcat(firstAddress,'comp',num2str(i),'.mat'));
        eng = load(strcat(firstAddress,'eng',num2str(i),'.mat'));
        
        F(i) = 0.33 * (minComp.minComp(i) ./ sum(computation.comp)) + 0.33 * (minEng.minEng(i) ./ sum(eng.eng)) + 0.33 * PDST.PDST(i);
    end
    
    save(strcat(firstAddress,'GF','.mat'),'F');
end
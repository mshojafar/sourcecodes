function [F] = goalFunction(firstAddress)

    computation = load(strcat(firstAddress,'comp','.mat'));
    eng = load(strcat(firstAddress,'eng','.mat'));
    PDST = load(strcat(firstAddress,'PDST','.mat'));
    minEng = load('MIterations\minEng.mat');
    minComp = load('MIterations\minComp.mat');
    
    
    F = 0.33 * (minComp.minComp ./ sum(computation.comp,2)) + 0.33 * (minEng.minEng ./ sum(eng.eng,2)) + 0.33 * PDST.PDST;

    save(strcat(firstAddress,'GF','.mat'),'F');
end
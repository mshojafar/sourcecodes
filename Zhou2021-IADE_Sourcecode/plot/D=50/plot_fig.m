clc;clear;
pool=[3,8,11,21];
for ii=1:4
    i=pool(ii);
    filepath1=strcat('./IADE/1_f',num2str(i));    
    datapath1=strcat(filepath1,'.mat');
    load(datapath1);
    plot(bestfitness,'LineWidth',1.5);
%    xticks(1e-15:1e5:1e14);
    xlabel('Evolution algebra(g)');
    ylabel('Function value');
%     legend('IADE');
    hold on;
    
    filepath2=strcat('./DE_rand_1/1_f',num2str(i));    
    datapath2=strcat(filepath2,'.mat');
    load(datapath2);
    plot(bestfitness,'LineWidth',1.5);
%     legend('DE/rand/1');
    hold on;
    
    filepath3=strcat('./DE_best_1/1_f',num2str(i));  
    datapath3=strcat(filepath3,'.mat');
    load(datapath3);
    plot(bestfitness,'LineWidth',1.5);
%     legend('DE/best/1');
    hold on;
    
    filepath4=strcat('./DE_rand2best_2/1_f',num2str(i));   
    datapath4=strcat(filepath4,'.mat');
    load(datapath4);
    plot(bestfitness,'LineWidth',1.5);
%     legend('DE/rand-to-best/2');
    hold on;
    
    legend('IADE','DE/rand/1','DE/best/1','DE/rand-to-best/2');

end
%IADE main program£¬

clc;clear;close all
global initial_flag
initial_flag = 0;
maxGen=5000;% Maximum iteration number
GEN=1;% Evolutionary algebra, or current iterative algebra

Dim=30;%Individual dimension, which can be modified to 10,30,50, 100
NP=100;%population size

Fmax=1;%scaling factor 
Fmin=0.4;

CRmin=0.6;%crossover rate 
CRmax=0.9;
fhd=str2func('cec17_func');

for iter=1:51  %Number of test
    for index=1:30%Index of test equation, different values correspond to different test functions
        t1 = [4 5 9 ];
        t2 = [1 2 3 8 10 11 12 13 14 15 16 17 18 20 21 22 23 24 25 26 27];
        t3 = [6 29 30];
        t4 = [7 19 28];
        
        if (ismember(index, t1))
            lb = -10;
            ub = 10;
        elseif (ismember(index, t2))
            lb = -100;
            ub = 100;
        elseif (ismember(index, t3))
            lb = -20;
            ub = 20;
        elseif (ismember(index, t4))
            lb = -50;
            ub = 50;
        end
        %% Initialize
        X=(ub-lb)*rand(NP,Dim)+lb; % Row X represents individual i, and column represents dimension j of individual i
        fitnessX=feval(fhd,X',index);
        [fitnessbestX,indexbestX]=min(fitnessX);
        bestX=X(indexbestX,:);%bestX means the position corresponding to the optimal value
        
        %% The iterative cycle
        for GEN=1:maxGen
            w=(cos((maxGen-GEN)/maxGen*pi)+1)/2;
            
            % Mutation
            F=Fmax-(Fmax-Fmin)*w;
            lamda=1-sqrt(GEN/maxGen);
            u=1-(GEN/maxGen)^2;
            V=mutation(X,bestX,F,lamda,u);
            %Boundary processing
            V=max(V,lb);
            V=min(V,ub);
            % Crossover
            CR=CRmin+(CRmax-CRmin)*w;
            [U,H]=crossover(X,V,CR);
            % Selection
            fitnessV=feval(fhd,V',index);
            fitnessU=feval(fhd,U',index);
            fitnessH=feval(fhd,H',index);
            
            for i=1:NP
                [~,out_x]=min([fitnessX(i);fitnessV(i);fitnessU(i);fitnessH(i)]);
                switch out_x
                    
                    case 1  %The target function for X is the smallest
                        Xavg=mean(X);
                        X(i,:)=(bestX+Xavg)./2;
                        
                    case 2  %The target function for V is the smallest
                        X(i,:)=V(i,:);
                        fitnessX(i)=fitnessV(i);
                        if fitnessV(i)<fitnessbestX
                            bestX=V(i,:);
                            fitnessbestX=fitnessV(i);
                        end
                        
                    case 3   %The target function for U is the smallest
                        X(i,:)=U(i,:);
                        fitnessX(i)=fitnessU(i);
                        if fitnessU(i)<fitnessbestX
                            bestX=U(i,:);
                            fitnessbestX=fitnessU(i);
                        end
                        
                    case 4   %The target function for H is the smallest
                        X(i,:)=H(i,:);
                        fitnessX(i)=fitnessH(i);
                        if fitnessH(i)<fitnessbestX
                            bestX=H(i,:);
                            fitnessbestX=fitnessH(i);
                        end
                    otherwise
                        error('ERROR!');
                end
            end
            % Recording
            fprintf('fun_mun %d-%d   %d   %f\n',index,iter,GEN,fitnessbestX);  %Output the optimal solution
            bestfitness(GEN)=fitnessbestX;
        end
        bestfit(index,iter)=fitnessbestX;%The final value of each run of each function
        groupFile=strcat('./D=30/',num2str(iter),'_f',num2str(index), '.mat');%Saving position
        save(groupFile,'bestfitness');
        save('./D=30/bestfit.mat','bestfit');
    end
end





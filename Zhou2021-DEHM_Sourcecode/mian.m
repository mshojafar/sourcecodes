%Main function£¬
clc;clear;close all

maxGen=5000;%Maximum iteration
GEN=1;%Evolutionary algebra, or current iteration algebra

Dim=30;%Individual dimensions, can be modified to 10,50,30
NP=100;%population size

Fmax=1;%scaling factor
Fmin=0.4; 

CRmin=0.6;%crossover rate
CRmax=0.9;
pool=[0.1*NP,0.7*NP,0.2*NP];
fhd=str2func('cec17_func');

for iter=1:51%Number of experiments
    for index=29:29
%Test equation index, different values correspond to different test functions
        
        [lb,ub,fobj] = Get_Functions_details(index);
        %% Initialization
        X=(ub-lb)*rand(NP,Dim)+lb; %  The X rows represent individual i, and the columns represent the dimension j of individual i
        fitnessX = (feval(fhd,X',index))';
        [fitnessbestX,indexbestX]=min(fitnessX);
        bestX=X(indexbestX,:);%bestX represents the location of the optimal value

        %% The iterative cycle
        for GEN=1:maxGen
            w=(cos((maxGen-GEN)/maxGen*pi)+1)/2;
            
            %Sorting
            fpop=sortrows([fitnessX,X]);
            X=fpop(:,2:Dim+1);
            
            %Orthogonal combination
            [finalx]=R_OED(X,fitnessX,bestX,index,fhd);
            
            %Quality population retention
            V=X;
            
            %General population variation operation (predicted position + Optimal position)
            temp1=0.2*NP+1;
            temp2=0.6*NP;
            flag=0;
            
            % variation
            F=Fmax-(Fmax-Fmin)*w;
            tempV=mutation(X,bestX,F,flag,finalx);%Call variable function
            V(temp1:temp2,:)=tempV(temp1:temp2,:);
            
            
            %Inferior population mutation operation (preferential learning)
            flag=1;
            temp1=0.8*NP+1;
            temp2=NP;
            tempV=mutation(X,bestX,F,flag,finalx);%Call variable function
            V(temp1:temp2,:)=tempV(temp1:temp2,:);
            
            %Boundary processing
            V=max(V,lb);
            V=min(V,ub);
            
            % Crossover
            CR=CRmin+(CRmax-CRmin)*w;
            [U]=crossover(X,V,CR);%Call crossover function
            
            % selection

            fitnessU=(feval(fhd,U',index))';
            fitnessV=(feval(fhd,V',index))';
            
            for i=1:NP
                [~,out_x]=min([fitnessX(i);fitnessV(i);fitnessU(i)]);
                
                switch out_x
                    
                    case 1  %The target function of X is the smallest
                        Xavg=mean(X);
                        X(i,:)=(bestX+Xavg)./2;
                        
                    case 2  %The target function of V is the smallest
                        X(i,:)=V(i,:);
                        fitnessX(i)=fitnessV(i);
                        if fitnessV(i)<fitnessbestX
                            bestX=V(i,:);
                            fitnessbestX=fitnessV(i);
                        end
                        
                    case 3   %The target function of U is the smallest
                        X(i,:)=U(i,:);
                        fitnessX(i)=fitnessU(i);
                        if fitnessU(i)<fitnessbestX
                            bestX=U(i,:);
                            fitnessbestX=fitnessU(i);
                        end
                        
                    otherwise
                        error('ERROR!');
                end
            end
            % Optimal value recording
            fprintf('fun_mun %d-%d   %d   %f\n',index,iter,GEN,fitnessbestX);  %Output the contemporary optimal solution
            bestfitness(GEN)=fitnessbestX;
        end
        bestfit(index,iter)=fitnessbestX;%The final value of each function for each run
        groupFile=strcat('./run',num2str(iter),'_f',num2str(index), '.mat');
        save(groupFile,'bestfitness');
        save('./bestfit.mat','bestfit');
    end
end




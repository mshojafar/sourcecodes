%%%%%%%%%%Esto se puede cambiar para que sea sensado cada 10 rounds o si no
%%%%%%%%%%que sea sensado siempre como se lo ha hecho en este caso si no
%%%%%%%%%%fuese asi se cambiaria  en la parte donde se elige los nodos CH y
%%%%%%%%%%se pondria q solo en tal rango sean sensados


clear;
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  %Field Dimensions - x and y maximum (in meters)
xm           =        100;
ym           =        100;
axis([0 xm 0 ym]);
grid on;
% for loop = 1:3
%     if loop==1
%         Eo           =        0.25;
%         x=0;
%         %         m=0.1;
%         %         a=1;
%     end
%     if loop==2
%         Eo           =        0.50;
%         x1=0;
%         %          m=0.1;
%         %         a=1;
%     end
%     if loop==3
%         Eo           =        1;
%         x2=0;
%         %       m=0.1;
%         %         a=1;
%     end
% 
% for loop = 1:3
%    if loop==1
%        m=0.1;
%        a=3;
%    end
%     if loop==2
%        m=0.2;
%        a=3;
%     end
%     if loop==3
%        m=0.3;
%         a=3;
%    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% VARIABLES  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n            =        100;          %Number of Nodes in the field
radio        =        n/10;         %Radio
m            =        0.1;          %Probbilidad de que un nodo se convierta en Cluster Head
p            =        0.1;          %Probabilidad
Eo           =        0.5;          %Initial Energy
ETX          =        50*(10^-9);   %Eelec=Etx=Erx
ERX          =        50*(10^-9);   %Eelec=Etx=Erx
Efs          =        10*(10^-12);  %Transmit Amplifier types
Emp          =        0.0013*(10^-12);
EDA          =        5*(10^-9);    %Data Aggregation Energy
a             =       1;            %\alpha - Percentage of nodes than are advanced - Values for Hetereogeneity
rmax         =        5000;         %maximum number of rounds
z            =        1;
do           =        sqrt(Efs/Emp);%distancia
rp           =        m*n+1;        %probabilidad de nodes advanced
j            =        0;            %variable para ordenar nodos advanced
u            =        0;
 x1=0;
 pv=0;
% %Thresholod for transmiting data to the cluster head
% t_h          =        10;           %%%%%%Hard Thres%%%%hold H(t)
% t_s          =        5;             %%%%%%Soft thres%%%%hold  S(t)
% sv           =        0;             %%%%%%previously Sensed value S(v)
% %methane level range
% l_r1         =        2.5;
% l_r          =        12.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%  END OF PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%  Creation of the random Sensor Network  %%%%%%%%%%%%%%
%figure(1);
for i=1:1:rp
    t=i^2;
    if (t+1)==(rp-1)
        rpr=(xm/(i));
        w=i;
        break;
    end
    if t>=(m*n+1)
        if ((m*n+1)-1)-t==0
            rpr=(xm/i);
            w=i;
        else
            rpr=(xm/(i));
            w=i;
        end
        
        break;
    end
end
for i=0:rpr:xm
    j=j+1;
    distancia1(j)=i;
end

l=j+1;
d=1;
rpr1=rpr;
sink.x=0.5*xm;
sink.y=0.5*ym;
S(n+1).xd=sink.x;
S(n+1).yd=sink.y;
for i=1:1:n
    
    S(i).xd=rand(1,1)*xm;
    XR(i)=S(i).xd;
    S(i).yd=rand(1,1)*ym;
    YR(i)=S(i).yd;
    S(i).G=0;
    ds(i)=sqrt((S(i).xd-sink.x)^2 + (S(i).yd-sink.y)^2 );
    
    
    % to avoid define normqal node much near to the sink so we re-define
    % normal node
    if ds(i)<10
    while  ds(i)<10
    S(i).xd=rand(1,1)*xm;
    XR(i)=S(i).xd;
    S(i).yd=rand(1,1)*ym;
    YR(i)=S(i).yd;
    S(i).G=0;
    ds(i)=sqrt((S(i).xd-S(n+1).xd)^2 + (S(i).yd-S(n+1).yd)^2 );       
    end
    end
    
    %initially there are no cluster heads only nodes
    S(i).type='N';
    temp_rnd0=i;
    
    %Random Election of Normal Nodes
    if (temp_rnd0>=m*n+1)
        S(i).E=Eo;
        S(i).ENERGY=0; % normal node indicator
        S(i).t=0;
        %%%%plot(S(i).xd,S(i).yd,'o');
        hold on;
    end
    
    %Random Election of Advanced Nodes
    if (temp_rnd0<m*n+1)
        S(i).E=Eo*(1+a);
        E1=S(i).E;
        S(i).ENERGY=1; % advanced node indicator
        S(i).t=0;
        %%%%plot(S(i).xd,S(i).yd,'+');
        if (mod(i,j)>=1) && (mod(i,j)<l)
            S(i).xd=rpr1;
            S(i).yd=distancia1(z);
            %%%%plot(S(i).xd,S(i).yd,'+');
            z=z+1;
        end
        % d part is used to place the last advanced node in the square
        % area.-
        if (mod(i,j)==0)
            S(i).xd=rpr1;
            S(i).yd=distancia1(z);
            z=1;
            d=d+1;
            rpr1=rpr*d;
        end
    end
end



%Primera Iteraccion
%figure(1);

%counter for CHs
countCHs=0;
%counter for CHs per round
rcountCHs=0;
cluster=1;
countCHs;
rcountCHs=rcountCHs+countCHs;
flag_first_dead=0;
 pn=0;
 pa=0;
for r=0:1:rmax
    r

%     if rem(r,3)==0
%       cv = l_r1 + (l_r-l_r1).*rand(1,1);  %Current sensing value C(v)
%       if (cv >5 && cv<10)  %periodo de normalidad
%             sv=0;
%       else if (cv>10)  %periodo de peligro
%                sv=1;
%           end
%       end  
%     end
  sum=0;
  % use this code to take a average of the energy of alive nodes in each
  % round
for j=1:1:n
if(S(j).E>0)
    sum=sum+S(j).E;
end
end

avg=sum/n;

    %Election Probability for Normal Nodes
    pnrm=( p/ (1+a*m));
    pn(r+1)=pnrm;
    %Election Probability for Advanced Nodes
    padv=( p*(1+a)/(1+a*m) )*(avg/E1);
    pa(r+1)=padv;
    
    %Operation for heterogeneous epoch
    % to avoid nodes two continues (adjacent rounds)times to be cluster head  
   if mod(r,2)==0
        if(mod(r, round(1/pnrm) )==0)
            for i=1:1:n
                S(i).G=0;
                S(i).cl=0;
            end
        end
    % end
   % if mod(r,2)==0
    %Operations for sub-epochs
    if(mod(r, round(1/padv) )==0)
        for i=1:1:n
            if(S(i).ENERGY==1)
                S(i).G=0;
                S(i).cl=0;
            end
        end
     end
   end
    
    hold off;
    
    %Numero de dead nodes
    dead=0;
    %Numero de dead Advanced Nodes
    dead_a=0;
    %Numero dead Normal Nodes
    dead_n=0;
    
    %Cuenta for bit transmitted to Bases Station and to Cluster Heads
    packets_TO_BS=0;
   % packets_TO_BS1=0;
    packets_TO_CH=0;
    %counter for bit transmitted to Bases Station and to Cluster Heads
    %per round
    PACKETS_TO_CH(r+1)=0;
    PACKETS_TO_BS(r+1)=0;
    %PACKETS_TO_BS1(r+1)=0;
    
 % figure(1);
    
    for i=1:1:n
        %checking if there is a dead node
        if (S(i).E<=0)
%         plot(S(i).xd,S(i).yd,'red .');
            dead=dead+1;
            if(S(i).ENERGY==1)   %%  si es =1 quiere decir q se mueve es un nodo adv
                dead_a=dead_a+1;
            end
            if(S(i).ENERGY==0) %%  si es =0 quiere decir q se no se mueve es un nodo normal
                dead_n=dead_n+1;
            end
            hold on;
        end
        % to check the node is alive
        if S(i).E>0
            S(i).type='N';   %%%% quiere decir q los nodos aun no estan muertos
            if (S(i).ENERGY==0)   % for checknig the normal node%%  dibuja los nodos normales con energia 0(identificacion)
 %          plot(S(i).xd,S(i).yd,'o','MarkerEdgeColor','b',...
  %                'MarkerFaceColor','w',...
%                   'MarkerSize',5);
            end
            if (S(i).ENERGY==1) % for checknig the advanced node
 %              plot(S(i).xd,S(i).yd,'D','MarkerEdgeColor','g',...
  %                 'MarkerFaceColor','g',...
%                   'MarkerSize',6);      
                %             circle(S(i).xd,S(i).yd,radio)
                hold on;
            end
        end
    end
 % plot(S(n+1).xd,S(n+1).yd,'-k*','MarkerSize',15,'MarkerEdgeColor','r');  %%  dibuja el nodo  sink
    
    
    STATISTICS(r+1).DEAD=dead;
    DEAD(r+1)=dead;
    DEAD_N(r+1)=dead_n;
    DEAD_A(r+1)=dead_a;
 
%cuando muere el primer nodo
    if (dead==1)
        if(flag_first_dead==0)
            first_dead=r;
            flag_first_dead=1;
        end
    end
    
countCHs=0;
cluster=1;  
sum1=0;
%  pa1=0;
%  pn1=0;
% STATISTICS(r+1).AVG=avg;
     pn1=(( pnrm / ( 1 - pnrm * mod(r,round(1/pnrm)))));
       pa1=(( padv / ( 1 - padv * mod(r,round(1/padv)))));
% INTIALIZATION FOR THE CLUSTER HEAD STRUCTUERING
for i=1:1:n
    
%       sum1=sum1+1;
%   sum=0;
for j=1:1:i
if(S(j).E>0)
    sum1=sum1+S(j).E;
end
end
% 
% avg=sum/i;

if (countCHs>=(m*n))
    break;
end
avg1=sum1/n;
% display(avg1);
% check the aliveness of the selected node
     if(S(i).E>0)
         
   temp_rand=rand; 
   
   % work on the alive node which doesn not have any set of the nodes which
   % are nominated as cluster head vicinity to it.
   
   if ( (S(i).G)<=0)
 
       
%        if  countCHs==0
%            
%            while (temp_rand<pnn1) || (temp_rand<paa1)              
%                paa1=pa1+temp_rand;
%                pnn1=pn1+temp_rand;
%                %                pn1=pn1+temp_rand;
%            end
%        end

       %%%%%%%%%%% THE MAIN IDEA OF THE THIS WORK%%%%%%%
       
       pn1=(( pnrm / ( 1 - pnrm * mod(r,round(1/pnrm))))*avg1);
       pa1=(( padv / ( 1 - padv * mod(r,round(1/padv))))*(1/avg1));
%         paa1=(( padv / ( 1 - padv * mod(r,round(1/padv)))));
           
             % check the node is the noraml noide and has probability to become a clustrer head   
            if    ( ( S(i).ENERGY==0 && ( temp_rand <= (((( pnrm / ( 1 - pnrm * mod(r,round(1/pnrm))))*avg1)))))) 
                %  if    ( ( S(i).ENERGY==0 && ( temp_rand <= (( pnrm / ( 1 - pnrm ))))))  ||    ( ( S(i).ENERGY==1 && ( temp_rand <= (( padv / ( 1 - padv ))))))
               
                % W is used to keep the coordinate of the cluster head in
                % each round to avoid emrging two neighbors nodes to be selected as cluster head  
                W(i).x=S(i).xd;
                W(i).y=S(i).yd;
                %                 if (mod(r,2)==0)
                % we sure that this node as a normal node can be selected as CH
                for j=1:1:i
                    distanciaz=sqrt( (S(i).xd-(W(j).x) )^2 + (S(i).yd-(W(j).y) )^2 );

                    if ((distanciaz>5))
                                                   % we sure that this node as a normal node can be selected as CH

                        countCHs=countCHs+1;
                        packets_TO_BS=packets_TO_BS+1;
                        PACKETS_TO_BS(r+1)=packets_TO_BS;
                        S(i).type='C';
                        S(i).G=100;
                        % save the cluster head coordinate in C
                        C(cluster).xd=S(i).xd;
                        C(cluster).yd=S(i).yd;
                        %  plot(S(i).xd,S(i).yd,'k*','MarkerEdgeColor','red',...
                        %    'MarkerFaceColor','red',...
                        %       'MarkerSize',6);
                        
                        % calculate the distance of the selected cluster
                        % head to the sink
                        distance=sqrt( (S(i).xd-(S(n+1).xd) )^2 + (S(i).yd-(S(n+1).yd) )^2 );
                        C(cluster).distance=distance;
                        C(cluster).id=i;
                        X(cluster)=S(i).xd;
                        Y(cluster)=S(i).yd;
                        cluster=cluster+1;
                        
                        %Calcula - Energy dissipated
                        distance;
                        %         if (sv==1)
                        if (distance>do)
                            S(i).E=S(i).E- ( (ETX+EDA)*(4000) + Emp*4000*( distance^4 ));
                        end
                        if (distance<=do)
                            S(i).E=S(i).E- ( (ETX+EDA)*(4000)  + Efs*4000*( distance^2 ));
                        end
                        %             packets_TO_BS1=packets_TO_BS1+1;
                        %             PACKETS_TO_BS1(r+1)=packets_TO_BS1;
                        
                        %         end
                        % each i we found our cluster head so no need to
                        % check again i° meighbors (means j) and we go to find next i
                        break;
                        
                        %                         end
                    end
                end
            end
%         end
           % check the node is the ADVANCED noide and has probability to become a clustrer head   

            if    ( ( S(i).ENERGY==1 && ( temp_rand <= (( (( padv / ( 1 - padv * mod(r,round(1/padv))))*(1/avg1)))))))
                %  if    ( ( S(i).ENERGY==0 && ( temp_rand <= (( pnrm / ( 1 - pnrm ))))))  ||    ( ( S(i).ENERGY==1 && ( temp_rand <= (( padv / ( 1 - padv ))))))
%                 W(i).x=S(i).xd;
%                 W(i).y=S(i).yd;
%                 %                 if (mod(r,2)==0)
%                 for j=1:1:i
%                     distanciaz=sqrt( (S(i).xd-(W(j).x) )^2 + (S(i).yd-(W(j).y) )^2 );
%                     if ((distanciaz>5))
                        countCHs=countCHs+1;
                        packets_TO_BS=packets_TO_BS+1;
                        PACKETS_TO_BS(r+1)=packets_TO_BS;
                        S(i).type='C';
                        S(i).G=100;
                        C(cluster).xd=S(i).xd;
                        C(cluster).yd=S(i).yd;
                        %  plot(S(i).xd,S(i).yd,'k*','MarkerEdgeColor','red',...
                        %    'MarkerFaceColor','red',...
                        %       'MarkerSize',6);
                        distance=sqrt( (S(i).xd-(S(n+1).xd) )^2 + (S(i).yd-(S(n+1).yd) )^2 );
                        C(cluster).distance=distance;
                        C(cluster).id=i;
                        X(cluster)=S(i).xd;
                        Y(cluster)=S(i).yd;
                        cluster=cluster+1;
                        
                        %Calcula - Energy dissipated
                        distance;
                        %         if (sv==1)
                        if (distance>do)
                            S(i).E=S(i).E- ( (ETX+EDA)*(4000) + Emp*4000*( distance^4 ));
                        end
                        if (distance<=do)
                            S(i).E=S(i).E- ( (ETX+EDA)*(4000)  + Efs*4000*( distance^2 ));
                        end
                        %             packets_TO_BS1=packets_TO_BS1+1;
                        %             PACKETS_TO_BS1(r+1)=packets_TO_BS1;
                        
                        %         end
%                         break;
                        
                        %                         end
%                     end
%                 end
            end
   end
     end
     
end


    % for showing the results in a plot
    STATISTICS(r+1).CLUSTERHEADS=cluster-1;
    CLUSTERHS(r+1)=cluster-1;
%     display(CLUSTERHS);
    % 
    % in the following 'for' we try to decrease the energy of the non-CH
    % nodes according to their distance of the nearest cluster head and
    % update their energy according to their distance to thier nearest
    % cluster heads
    %Cuantos se convierten en Cluster Head en cada round  (CLUSTERHS)
    %Election of Associated Cluster Head for Normal Nodes  Numero de normal
    %nodes que tiene a cargo el cluster head
    for i=1:1:n
    if ( S(i).type=='N' && S(i).E>0 )
        if(cluster-1>=1)
            min_dis=sqrt( (S(i).xd-S(n+1).xd)^2 + (S(i).yd-S(n+1).yd)^2 );
            min_dis_cluster=1;
            for c=1:1:cluster-1
                temp=min(min_dis,sqrt( (S(i).xd-C(c).xd)^2 + (S(i).yd-C(c).yd)^2 ) );
                if ( temp<min_dis )
                    min_dis=temp;
                    min_dis_cluster=c;
                end
            end
                
          %Energy dissipated by associated Cluster Head
            min_dis;
            if (min_dis>do)
                S(i).E=S(i).E- ( ETX*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis));
            end
            if (min_dis<=do)
                S(i).E=S(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis));
            end
            %Energy dissipated
            if(min_dis>0)
                S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E- ( (ERX + EDA)*4000 );
                PACKETS_TO_CH(r+1)=n-dead-cluster+1;
            end
            
            S(i).min_dis=min_dis;
            S(i).min_dis_cluster=min_dis_cluster;
            
        end
    end
end
    if DEAD(r+1)==n
        rmax=r+1;
        break;
    end
    
    hold on;
    
    countCHs;
  
    rcountCHs=rcountCHs+countCHs;
    sum=0;
for i=1:1:n
if(S(i).E>0)
    sum=sum+S(i).E;
end
end

avg=sum/n;
STATISTICS(r+1).AVG=avg;
uu(r+1).AVG=(STATISTICS(r+1).AVG)*100;
sum;
    grid on;
%     if Eo==0.25
%         x(r+1)=uu(r+1).AVG;
%     end
  
        x1(r+1)=uu(r+1).AVG;

%     if Eo==1
%         x2(r+1)=uu(r+1).AVG;
%     end
%     [vx,vy]=voronoi(X,Y);
% plot(X,Y,'r*',vx,vy,'b-');
% hold on;
% voronoi(X,Y);
% axis([0 xm 0 ym]);
end % round


alive=u;
if m==0.1
    x=u;
end
% if m==0.2
%     x1=u;
% end
% if m==0.3
%     x2=u;
% end
for i=1:1:rmax
    alive(i)=n-DEAD(i);
%     if m==0.1
%         x(i)= alive(i);
%     end
% if m==0.2
%     x1(i)= alive(i);
% end
% if m==0.3
%     x2(i)= alive(i);
% end 
   
end
pack=PACKETS_TO_BS;
save  pack.mat;
save  alive.mat;
save  x1.mat;
save  pa.mat;
save  pn.mat;

% end
% figure(3)
% r=1:1:rmax;
% plot(x(r),r,'--',x1(r),r,'g',x2(r),r,':');
% legend(['m=0.2','  ','\alpha=3','  ','Eo=0.25'],['m=0.2','  ','\alpha=3','  ','Eo=0.50'],['m=0.2','  ','\alpha=3','  ','Eo=1']);
%     xlabel('Average Energy of Each Node');
%     ylabel('Round Number');   
%      hold on;

% figure(4)
% for 
% end
% figure(3)
% r=0:1:rmax;
%     plot(x,r,'--rs',x1,r,'g',x2,r,':');
%     legend(['m=0.1','  ','\alpha=1','  ','Eo=0.25'],['m=0.1','  ','\alpha=1','  ','Eo=0.50'],['m=0.1','  ','\alpha=1','  ','Eo=1']);
%     xlabel('Average Energy of Each Node');
%     ylabel('Round Number');
    

% figure(3)
% r=0:1:rmax;
% plot(x,r,'--rs',x1,r,'g',x2,r,':');
% xlabel('Average Energy of Each Node');
% ylabel('Round Number');
% legend(['m=0.1','  ','\alpha=1','  ','Eo=0.25'],['m=0.1','  ','\alpha=1','  ','Eo=0.50'],['m=0.1','  ','\alpha=1','  ','Eo=1']);
%      hold on;

% figure(3)
% r=0:1:rmax;
%        plot(x,r,'--rs',x1,r,'g',x2,r,':');
%     xlabel('Average Energy of Each Node');
%     ylabel('Round Number');
% %   legend(['Eo=' num2str(Eo),'','Eo=' num2str(Eo),'','Eo=' num2str(Eo)]);
%      hold on;
 
 
% 
% figure(4)
% r=1:1:rmax;
% plot(r,alive,'--');
% legend(['m=0.1','  ','\alpha=1']);
% ylabel('Number of AliveNodes')
% xlabel('Round Number')
% figure(4)
% r=1:1:rmax;
% plot(r,x,'--',r,x1,'g',r,x2,':');
% legend(['m=0.1','  ','\alpha=3'],['m=0.2','  ','\alpha=3'],['m=0.3','  ','\alpha=3']);
% ylabel('Number of AliveNodes')
% xlabel('Round Number')

% legend('SEP','Location','northeast');
% end
% end
 %  figure(4)
%  for r=0:1:rmax;
%     plot([r r+1],[STATISTICS(r+1).AVG STATISTICS(r+2).AVG],'red');
%     plot(r,alive,'color','blue');
%     hold on;
%  end


%%
%%SEP ORIGINAL

clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Field Dimensions - x and y maximum (in meters)
xm=100;
ym=100;

%x and y Coordinates of the Sink
sink.x=0.5*xm;
sink.y=0.5*ym;

%Number of Nodes in the field
n=100;
u=0;
%Optimal Election Probability of a node
%to become cluster head
p=0.1;

%Energy Model (all values in Joules)
%Initial Energy 
Eo=0.5;
%Eelec=Etx=Erx
ETX=50*0.000000001;
ERX=50*0.000000001;
%Transmit Amplifier types
Efs=10*0.000000000001;
Emp=0.0013*0.000000000001;
%Data Aggregation Energy
EDA=5*0.000000001;

%Values for Hetereogeneity
%Percentage of nodes than are advanced
m=0.1;
%\alpha
a=1;

%maximum number of rounds
rmax=5000;

%%%%%%%%%%%%%%%%%%%%%%%%% END OF PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%

%Computation of do
do=sqrt(Efs/Emp); 

%Creation of the random Sensor Network
%este figure(1);
for i=1:1:n
    S(i).xd=rand(1,1)*xm;
    XR(i)=S(i).xd;
    S(i).yd=rand(1,1)*ym;
    YR(i)=S(i).yd;
    S(i).G=0;
    %initially there are no cluster heads only nodes
    S(i).type='N';
   
    temp_rnd0=i;
    %Random Election of Normal Nodes
    if (temp_rnd0>=m*n+1) 
        S(i).E=Eo;
        S(i).ENERGY=0;
        %%%%plot(S(i).xd,S(i).yd,'o');
        hold on;
    end
    %Random Election of Advanced Nodes
    if (temp_rnd0<m*n+1)  
        S(i).E=Eo*(1+a);
        S(i).ENERGY=1;
        %%%%plot(S(i).xd,S(i).yd,'+');
         hold on;
    end
end

S(n+1).xd=sink.x;
S(n+1).yd=sink.y;
%%%%plot(S(n+1).xd,S(n+1).yd,'x');
    
        
%First Iteration
%estefigure(1);

%counter for CHs
countCHs=0;
%counter for CHs per round
rcountCHs=0;
cluster=1;

countCHs;
rcountCHs=rcountCHs+countCHs;
flag_first_dead=0;

for r=0:1:rmax
   % r

  %Election Probability for Normal Nodes
  pnrm=( p/ (1+a*m) );
  %Election Probability for Advanced Nodes
  padv= ( p*(1+a)/(1+a*m) );
    
  %Operation for heterogeneous epoch
  if(mod(r, round(1/pnrm) )==0)
    for i=1:1:n
        S(i).G=0;
        S(i).cl=0;
    end
  end

 %Operations for sub-epochs
 if(mod(r, round(1/padv) )==0)
    for i=1:1:n
        if(S(i).ENERGY==1)
            S(i).G=0;
            S(i).cl=0;
        end
    end
  end

 
hold off;

%Number of dead nodes
dead=0;
%Number of dead Advanced Nodes
dead_a=0;
%Number of dead Normal Nodes
dead_n=0;

%counter for bit transmitted to Bases Station and to Cluster Heads
packets_TO_BS=0;
packets_TO_CH=0;
%counter for bit transmitted to Bases Station and to Cluster Heads 
%per round
PACKETS_TO_CH(r+1)=0;
PACKETS_TO_BS(r+1)=0;

%estefigure(1);

for i=1:1:n
    %checking if there is a dead node
    if (S(i).E<=0)
 %este       plot(S(i).xd,S(i).yd,'red .');
        dead=dead+1;
        if(S(i).ENERGY==1)
            dead_a=dead_a+1;
        end
        if(S(i).ENERGY==0)
            dead_n=dead_n+1;
        end
        hold on;    
    end
    if S(i).E>0
        S(i).type='N';
        if (S(i).ENERGY==0)  
 %este       plot(S(i).xd,S(i).yd,'o');
        end
        if (S(i).ENERGY==1)  
  %este      plot(S(i).xd,S(i).yd,'+');
        end
        hold on;
    end
end
%este plot(S(n+1).xd,S(n+1).yd,'x');


STATISTICS(r+1).DEAD=dead;
DEAD(r+1)=dead;
DEAD_N(r+1)=dead_n;
DEAD_A(r+1)=dead_a;

%When the first node dies
if (dead==1)
    if(flag_first_dead==0)
        first_dead=r;
        flag_first_dead=1;
    end
end

countCHs=0;
cluster=1;
for i=1:1:n
   if(S(i).E>0)
   temp_rand=rand;     
   if ( (S(i).G)<=0)

 %Election of Cluster Heads for normal nodes
 if( ( S(i).ENERGY==0 && ( temp_rand <= ( pnrm / ( 1 - pnrm * mod(r,round(1/pnrm)) )) ) )  )
     
     countCHs=countCHs+1;
     packets_TO_BS=packets_TO_BS+1;
     PACKETS_TO_BS(r+1)=packets_TO_BS;
     
     S(i).type='C';
     S(i).G=100;
     C(cluster).xd=S(i).xd;
     C(cluster).yd=S(i).yd;
 %este    plot(S(i).xd,S(i).yd,'k*');
     
     distance=sqrt( (S(i).xd-(S(n+1).xd) )^2 + (S(i).yd-(S(n+1).yd) )^2 );
     C(cluster).distance=distance;
     C(cluster).id=i;
     X(cluster)=S(i).xd;
     Y(cluster)=S(i).yd;
     cluster=cluster+1;
     
     %Calculation of Energy dissipated
     distance;
     if (distance>do)
         S(i).E=S(i).E- ( (ETX+EDA)*(4000) + Emp*4000*( distance*distance*distance*distance ));
     end
     if (distance<=do)
         S(i).E=S(i).E- ( (ETX+EDA)*(4000)  + Efs*4000*( distance * distance ));
     end
 end
 


 %Election of Cluster Heads for Advanced nodes
 if( ( S(i).ENERGY==1 && ( temp_rand <= ( padv / ( 1 - padv * mod(r,round(1/padv)) )) ) )  )
        
            countCHs=countCHs+1;
            packets_TO_BS=packets_TO_BS+1;
            PACKETS_TO_BS(r+1)=packets_TO_BS;
            
            S(i).type='C';
            S(i).G=100;
            C(cluster).xd=S(i).xd;
            C(cluster).yd=S(i).yd;
 %este           plot(S(i).xd,S(i).yd,'k*');
            
            distance=sqrt( (S(i).xd-(S(n+1).xd) )^2 + (S(i).yd-(S(n+1).yd) )^2 );
            C(cluster).distance=distance;
            C(cluster).id=i;
            X(cluster)=S(i).xd;
            Y(cluster)=S(i).yd;
            cluster=cluster+1;
            
            %Calculation of Energy dissipated
            distance;
            if (distance>do)
                S(i).E=S(i).E- ( (ETX+EDA)*(4000) + Emp*4000*( distance*distance*distance*distance )); 
            end
            if (distance<=do)
                S(i).E=S(i).E- ( (ETX+EDA)*(4000)  + Efs*4000*( distance * distance )); 
            end
        end     
    
    end
  end 
end



STATISTICS(r+1).CLUSTERHEADS=cluster-1;
CLUSTERHS(r+1)=cluster-1;

%Election of Associated Cluster Head for Normal Nodes
for i=1:1:n
    if ( S(i).type=='N' && S(i).E>0 )
        if(cluster-1>=1)
            min_dis=sqrt( (S(i).xd-S(n+1).xd)^2 + (S(i).yd-S(n+1).yd)^2 );
            min_dis_cluster=1;
            for c=1:1:cluster-1
                temp=min(min_dis,sqrt( (S(i).xd-C(c).xd)^2 + (S(i).yd-C(c).yd)^2 ) );
                if ( temp<min_dis )
                    min_dis=temp;
                    min_dis_cluster=c;
                end
            end
            
            %Energy dissipated by associated Cluster Head
            min_dis;
            if (min_dis>do)
                S(i).E=S(i).E- ( ETX*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis));
            end
            if (min_dis<=do)
                S(i).E=S(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis));
            end
            %Energy dissipated
            if(min_dis>0)
                S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E- ( (ERX + EDA)*4000 );
                PACKETS_TO_CH(r+1)=n-dead-cluster+1;
            end
            
            S(i).min_dis=min_dis;
            S(i).min_dis_cluster=min_dis_cluster;
            
        end
    end
end
hold on;
%     if DEAD(r+1)==n
%         rmax=r+1;
%         break;
%     end
countCHs;
rcountCHs=rcountCHs+countCHs;
  sum=0;
for i=1:1:n
if(S(i).E>0)
    sum=sum+S(i).E;
end
end
avg=sum/n;
STATISTICS(r+1).AVG=avg;
uu(r+1).AVG=(STATISTICS(r+1).AVG)*100;
sum;
alive1=u;
x2(r+1)=uu(r+1).AVG;

% if m==0.1
%     x=u;
% end
% if m==0.2
%     x1=u;
% end
% if m==0.3
%     x2=u;
% end

%Code for Voronoi Cells
%Unfortynately if there is a small
%number of cells, Matlab's voronoi
%procedure has some problems

%[vx,vy]=voronoi(X,Y);
%plot(X,Y,'r*',vx,vy,'b-');
% hold on;
% voronoi(X,Y);
% axis([0 xm 0 ym]);
   

end
for i=1:1:rmax
    alive1(i)=n-DEAD(i);
%     if m==0.1
%         x(i)= alive(i);
%     end
% if m==0.2
%     x1(i)= alive(i);
% end
% if m==0.3
%     x2(i)= alive(i);
% end 
   
end
pack1=PACKETS_TO_BS;
save  pack1.mat;
save alive1.mat;
save  x2.mat;
% figure(3)
% r=1:1:rmax;
% plot(r,alive,'--');
% legend(['m=0.1','  ','\alpha=1']);
% ylabel('Number of AliveNodes')
% xlabel('Round Number')
% hold on;
%%
%%LEACH-DCHS

clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Field Dimensions - x and y maximum (in meters)
xm=100;
ym=100;

%x and y Coordinates of the Sink
sink.x=0.5*xm;
sink.y=0.5*ym;

%Number of Nodes in the field
n=100;

%Optimal Election Probability of a node
%to become cluster head
p=0.1;

%Energy Model (all values in Joules)
%Initial Energy 
Eo=0.5;
%Eelec=Etx=Erx
ETX=50*0.000000001;
ERX=50*0.000000001;
%Transmit Amplifier types
Efs=10*0.000000000001;
Emp=0.0013*0.000000000001;
%Data Aggregation Energy
EDA=5*0.000000001;

%Values for Hetereogeneity
%Percentage of nodes than are advanced
m=0.1;
%\alpha
a=1;
u=0;
x3=0;
r_s=0;
%maximum number of rounds
rmax=5000;

%%%%%%%%%%%%%%%%%%%%%%%%% END OF PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%

%Computation of do
do=sqrt(Efs/Emp);

%Creation of the random Sensor Network
figure(1);
for i=1:1:n
    S(i).xd=rand(1,1)*xm;
    XR(i)=S(i).xd;
    S(i).yd=rand(1,1)*ym;
    YR(i)=S(i).yd;
    S(i).G=0;
    %initially there are no cluster heads only nodes
    S(i).type='N';
   
    temp_rnd0=i;
    %Random Election of Normal Nodes
    if (temp_rnd0>=m*n+1) 
        S(i).E=Eo;
        S(i).ENERGY=0;
  %      plot(S(i).xd,S(i).yd,'o');
        hold on;
    end
    %Random Election of Advanced Nodes
    if (temp_rnd0<m*n+1)  
        S(i).E=Eo*(1+a);
        S(i).ENERGY=1;
   %     plot(S(i).xd,S(i).yd,'+');
        hold on;
    end
end

S(n+1).xd=sink.x;
S(n+1).yd=sink.y;
%plot(S(n+1).xd,S(n+1).yd,'x');
    
        
%First Iteration
figure(1);

%counter for CHs
countCHs=0;
%counter for CHs per round
rcountCHs=0;
cluster=1;

countCHs;
rcountCHs=rcountCHs+countCHs;
flag_first_dead=0;

for r=0:1:rmax
    r

  %Operation for epoch
  if(mod(r, round(1/p) )==0)
    for i=1:1:n
        S(i).G=0;
        S(i).cl=0;
    end
  end

hold off;

%Number of dead nodes
dead=0;
%Number of dead Advanced Nodes
dead_a=0;
%Number of dead Normal Nodes
dead_n=0;

%counter for bit transmitted to Bases Station and to Cluster Heads
packets_TO_BS=0;
packets_TO_CH=0;
%counter for bit transmitted to Bases Station and to Cluster Heads 
%per round
PACKETS_TO_CH(r+1)=0;
PACKETS_TO_BS(r+1)=0;

figure(1);

for i=1:1:n
    %checking if there is a dead node
    if (S(i).E<=0)
 %       plot(S(i).xd,S(i).yd,'red .');
        dead=dead+1;
        if(S(i).ENERGY==1)
            dead_a=dead_a+1;
        end
        if(S(i).ENERGY==0)
            dead_n=dead_n+1;
        end
        hold on;    
    end
    if S(i).E>0
        S(i).type='N';
        if (S(i).ENERGY==0)  
 %       plot(S(i).xd,S(i).yd,'o');
        end
        if (S(i).ENERGY==1)  
%        plot(S(i).xd,S(i).yd,'+');
        end
        hold on;
    end
end
%plot(S(n+1).xd,S(n+1).yd,'x');


STATISTICS(r+1).DEAD=dead;
DEAD(r+1)=dead;
DEAD_N(r+1)=dead_n;
DEAD_A(r+1)=dead_a;

%When the first node dies
if (dead==1)
    if(flag_first_dead==0)
        first_dead=r;
        flag_first_dead=1;
    end
end

countCHs=0;
cluster=1;
for i=1:1:n
  r_s=r_s+1;
   if(S(i).E>0)
   temp_rand=rand;     
   if ( (S(i).G)<=0)

 %Election of Cluster-Heads
 if(temp_rand<= (p/(1-p*mod(r,round(1/p))))*(((S(i).E/Eo)+(r_s/(1/p))*(1-(S(i).E/Eo)))))
     %*((S(i).E/Eo)+(r_s/(1/p))*(1-(S(i).E/Eo)))
     countCHs=countCHs+1;
     packets_TO_BS=packets_TO_BS+1;
     PACKETS_TO_BS(r+1)=packets_TO_BS;
     
     S(i).type='C';
     S(i).G=round(1/p)-1;
     C(cluster).xd=S(i).xd;
     C(cluster).yd=S(i).yd;
%     plot(S(i).xd,S(i).yd,'k*');
     
     distance=sqrt( (S(i).xd-(S(n+1).xd) )^2 + (S(i).yd-(S(n+1).yd) )^2 );
     C(cluster).distance=distance;
     C(cluster).id=i;
     X(cluster)=S(i).xd;
     Y(cluster)=S(i).yd;
     cluster=cluster+1;
     
     %Calculation of Energy dissipated
     distance;
     if (distance>do)
         S(i).E=S(i).E- ( (ETX+EDA)*(4000) + Emp*4000*( distance*distance*distance*distance ));
     end
     if (distance<=do)
         S(i).E=S(i).E- ( (ETX+EDA)*(4000)  + Efs*4000*( distance * distance ));
     end
 end
    
    end
   end 
 r_s=0;
end

STATISTICS(r+1).CLUSTERHEADS=cluster-1;
CLUSTERHS(r+1)=cluster-1;

%Election of Associated Cluster Head for Normal Nodes
for i=1:1:n
   if ( S(i).type=='N' && S(i).E>0 )
     if(cluster-1>=1)
       min_dis=sqrt( (S(i).xd-S(n+1).xd)^2 + (S(i).yd-S(n+1).yd)^2 );
       min_dis_cluster=1;
       for c=1:1:cluster-1
           temp=min(min_dis,sqrt( (S(i).xd-C(c).xd)^2 + (S(i).yd-C(c).yd)^2 ) );
           if ( temp<min_dis )
               min_dis=temp;
               min_dis_cluster=c;
           end
       end
       
       %Energy dissipated by associated Cluster Head
            min_dis;
            if (min_dis>do)
                S(i).E=S(i).E- ( ETX*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis)); 
            end
            if (min_dis<=do)
                S(i).E=S(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis)); 
            end
        %Energy dissipated
        if(min_dis>0)
          S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E- ( (ERX + EDA)*4000 ); 
         PACKETS_TO_CH(r+1)=n-dead-cluster+1; 
        end

       S(i).min_dis=min_dis;
       S(i).min_dis_cluster=min_dis_cluster;
           
   end
 end
end
hold on;
% if DEAD(r+1)==n
%     rmax=r+1;
%     break;
% end
countCHs;
rcountCHs=rcountCHs+countCHs;
  sum=0;
for i=1:1:n
if(S(i).E>0)
    sum=sum+S(i).E;
end
end

avg=sum/n;
STATISTICS(r+1).AVG=avg;
uu(r+1).AVG=(STATISTICS(r+1).AVG)*100;
sum;
x3(r+1)=uu(r+1).AVG;


%Code for Voronoi Cells
%Unfortynately if there is a small
%number of cells, Matlab's voronoi
%procedure has some problems

%[vx,vy]=voronoi(X,Y);
%plot(X,Y,'r*',vx,vy,'b-');
% hold on;
% voronoi(X,Y);
% axis([0 xm 0 ym]);

end

alive2=u;

for i=1:1:rmax
    alive2(i)=n-DEAD(i);
end
 pack2=PACKETS_TO_BS;
save  pack2.mat;
save alive2.mat;
save  x3.mat;
% figure(5)
% r= 1:1:rmax;
% plot(r,alive,'r');  



%%

% %%A-LEACH
% 
% clear all;
% clc;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %Field Dimensions - x and y maximum (in meters)
% xm=100;
% ym=100;
% 
% %x and y Coordinates of the Sink
% sink.x=0.5*xm;
% sink.y=0.5*ym;
% 
% %Number of Nodes in the field
% n=100;
% 
% %Optimal Election Probability of a node
% %to become cluster head
% p=0.1;
% 
% %Energy Model (all values in Joules)
% %Initial Energy 
% Eo=0.5;
% %Eelec=Etx=Erx
% ETX=50*0.000000001;
% ERX=50*0.000000001;
% %Transmit Amplifier types
% Efs=10*0.000000000001;
% Emp=0.0013*0.000000000001;
% %Data Aggregation Energy
% EDA=5*0.000000001;
% 
% %Values for Hetereogeneity
% %Percentage of nodes than are advanced
% m=0.1;
% %\alpha
% a=1;
% u=0;
% x4=0;
% %maximum number of rounds
% rmax=5000;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%% END OF PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%
% 
% %Computation of do
% do=sqrt(Efs/Emp);
% 
% %Creation of the random Sensor Network
% %figure(1);
% for i=1:1:n
%     S(i).xd=rand(1,1)*xm;
%     XR(i)=S(i).xd;
%     S(i).yd=rand(1,1)*ym;
%     YR(i)=S(i).yd;
%     S(i).G=0;
%     %initially there are no cluster heads only nodes
%     S(i).type='N';
%    
%     temp_rnd0=i;
%     %Random Election of Normal Nodes
%     if (temp_rnd0>=m*n+1) 
%         S(i).E=Eo;
%         S(i).ENERGY=0;
%   %      plot(S(i).xd,S(i).yd,'o');
%         hold on;
%     end
%     %Random Election of Advanced Nodes
%     if (temp_rnd0<m*n+1)  
%         S(i).E=Eo*(1+a);
%         S(i).ENERGY=1;
%    %     plot(S(i).xd,S(i).yd,'+');
%         hold on;
%     end
% end
% 
% S(n+1).xd=sink.x;
% S(n+1).yd=sink.y;
% %plot(S(n+1).xd,S(n+1).yd,'x');
%     
%         
% %First Iteration
% %figure(1);
% 
% %counter for CHs
% countCHs=0;
% %counter for CHs per round
% rcountCHs=0;
% cluster=1;
% 
% countCHs;
% rcountCHs=rcountCHs+countCHs;
% flag_first_dead=0;
% 
% for r=0:1:rmax
%     r
% 
%   %Operation for epoch
%   if(mod(r, round(1/p) )==0)
%     for i=1:1:n
%         S(i).G=0;
%         S(i).cl=0;
%     end
%   end
% 
% hold off;
% 
% %Number of dead nodes
% dead=0;
% %Number of dead Advanced Nodes
% dead_a=0;
% %Number of dead Normal Nodes
% dead_n=0;
% 
% %counter for bit transmitted to Bases Station and to Cluster Heads
% packets_TO_BS=0;
% packets_TO_CH=0;
% %counter for bit transmitted to Bases Station and to Cluster Heads 
% %per round
% PACKETS_TO_CH(r+1)=0;
% PACKETS_TO_BS(r+1)=0;
% 
% %figure(1);
% 
% for i=1:1:n
%     %checking if there is a dead node
%     if (S(i).E<=0)
%  %      plot(S(i).xd,S(i).yd,'red .');
%         dead=dead+1;
%         if(S(i).ENERGY==1)
%             dead_a=dead_a+1;
%         end
%         if(S(i).ENERGY==0)
%             dead_n=dead_n+1;
%         end
%         hold on;    
%     end
%     if S(i).E>0
%         S(i).type='N';
%         if (S(i).ENERGY==0)  
%   %     plot(S(i).xd,S(i).yd,'o');
%         end
%         if (S(i).ENERGY==1)  
%    %     plot(S(i).xd,S(i).yd,'+');
%         end
%         hold on;
%     end
% end
% %plot(S(n+1).xd,S(n+1).yd,'x');
% 
% 
% STATISTICS(r+1).DEAD=dead;
% DEAD(r+1)=dead;
% DEAD_N(r+1)=dead_n;
% DEAD_A(r+1)=dead_a;
% 
% %When the first node dies
% if (dead==1)
%     if(flag_first_dead==0)
%         first_dead=r;
%         flag_first_dead=1;
%     end
% end
% 
% countCHs=0;
% cluster=1;
%  %*r_s=0;
%  cl1=0;
% for i=1:1:n 
%    %*  r_s=r_s+1;
%    if(S(i).E>0)
%    temp_rand=rand;     
%    if ( (S(i).G)<=0)
% 
%  %Election of Cluster-Heads
%  if(temp_rand<= (p/(n-p*mod(r,round(n/p))))+((S(i).E/Eo)*(p/n)))    
%    %*((S(i).E/Eo)+(r_s/(1/p))*(1-(S(i).E/Eo)))
%    cl1=cl1+1;
%      countCHs=countCHs+1;
%      packets_TO_BS=packets_TO_BS+1;
%      PACKETS_TO_BS(r+1)=packets_TO_BS;
%      
%      S(i).type='C';
%      S(i).G=round(1/p)-1;
%      C(cluster).xd=S(i).xd;
%      C(cluster).yd=S(i).yd;
%   %   plot(S(i).xd,S(i).yd,'k*');
%      
%      distance=sqrt( (S(i).xd-(S(n+1).xd) )^2 + (S(i).yd-(S(n+1).yd) )^2 );
%      C(cluster).distance=distance;
%      C(cluster).id=i;
%      X(cluster)=S(i).xd;
%      Y(cluster)=S(i).yd;
%      cluster=cluster+1;
%      
%      %Calculation of Energy dissipated
%      distance;
%      if (distance>do)
%          S(i).E=S(i).E- ( (ETX+EDA)*(4000) + Emp*4000*( distance*distance*distance*distance ));
%      end
%      if (distance<=do)
%          S(i).E=S(i).E- ( (ETX+EDA)*(4000)  + Efs*4000*( distance * distance ));
%      end
%  end
%     
%     end
%    end 
%  % r_s=0;
% 
% end
% 
% STATISTICS(r+1).CLUSTERHEADS=cluster-1;
% CLUSTERHS(r+1)=cluster-1;
% 
% %Election of Associated Cluster Head for Normal Nodes
% for i=1:1:n
%    if ( S(i).type=='N' && S(i).E>0 )
%      if(cluster-1>=1)
%        min_dis=sqrt( (S(i).xd-S(n+1).xd)^2 + (S(i).yd-S(n+1).yd)^2 );
%        min_dis_cluster=1;
%        for c=1:1:cluster-1
%            temp=min(min_dis,sqrt( (S(i).xd-C(c).xd)^2 + (S(i).yd-C(c).yd)^2 ) );
%            if ( temp<min_dis )
%                min_dis=temp;
%                min_dis_cluster=c;
%            end
%        end
%        
%        %Energy dissipated by associated Cluster Head
%             min_dis;
%             if (min_dis>do)
%                 S(i).E=S(i).E- ( ETX*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis)); 
%             end
%             if (min_dis<=do)
%                 S(i).E=S(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis)); 
%             end
%         %Energy dissipated
%         if(min_dis>0)
%           S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E- ( (ERX + EDA)*4000 ); 
%          PACKETS_TO_CH(r+1)=n-dead-cluster+1; 
%         end
% 
%        S(i).min_dis=min_dis;
%        S(i).min_dis_cluster=min_dis_cluster;
%            
%    end
%  end
% end
% hold on;
% % if DEAD(r+1)==n
% %     rmax=r+1;
% %     break;
% % end
% countCHs;
% rcountCHs=rcountCHs+countCHs;
%   sum=0;
% for i=1:1:n
% if(S(i).E>0)
%     sum=sum+S(i).E;
% end
% end
% 
% 
% 
% 
% %Code for Voronoi Cells
% %Unfortynately if there is a small
% %number of cells, Matlab's voronoi
% %procedure has some problems
% 
% %[vx,vy]=voronoi(X,Y);
% %plot(X,Y,'r*',vx,vy,'b-');
% % hold on;
% % voronoi(X,Y);
% % axis([0 xm 0 ym]);
% 
% end
% pack3=PACKETS_TO_BS;
% alive3=u;
% 
% 
% j=0;
% CL=0;
% for i=1:rmax
%     if CLUSTERHS(i)~=0
%         j=j+1;
%         CL(j)=CLUSTERHS(i);  
%         DEAD1(j)=DEAD(i);
%         avg=sum/n;
% STATISTICS(i+1).AVG=avg;
% uu(i+1).AVG=(STATISTICS(i+1).AVG)*100;
% sum;
% x4(i+1)=uu(i+1).AVG;
%     end
%     
% end
% 
% for i=1:1:j
%     alive3(i)=n-DEAD1(i);
% end
% 
% for i=(j+1):1:rmax
%     alive3(i)=0;
%     x4(i)=0;
%     pack3(i)=0;
% end
% 
% save  pack3.mat;
% save alive3.mat;
% save  x4.mat;
% % figure(5)
% % r= 1:1:rmax;
% % plot(r,alive,'r');  
%%
%%M-SEP ORIGINAL

clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Field Dimensions - x and y maximum (in meters)
xm=100;
ym=100;

%x and y Coordinates of the Sink
sink.x=0.5*xm;
sink.y=0.5*ym;

%Number of Nodes in the field
n=100;
u=0;
%Optimal Election Probability of a node
%to become cluster head
p=0.1;

%Energy Model (all values in Joules)
%Initial Energy 
Eo=0.5;
%Eelec=Etx=Erx
ETX=50*0.000000001;
ERX=50*0.000000001;
%Transmit Amplifier types
Efs=10*0.000000000001;
Emp=0.0013*0.000000000001;
%Data Aggregation Energy
EDA=5*0.000000001;

%Values for Hetereogeneity
%Percentage of nodes than are advanced
m=0.1;
%\alpha
a=1;

%maximum number of rounds
rmax=5000;

%%%%%%%%%%%%%%%%%%%%%%%%% END OF PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%

%Computation of do
do=sqrt(Efs/Emp); 

%Creation of the random Sensor Network
%este figure(1);
for i=1:1:n
    S(i).xd=rand(1,1)*xm;
    XR(i)=S(i).xd;
    S(i).yd=rand(1,1)*ym;
    YR(i)=S(i).yd;
    S(i).G=0;
    %initially there are no cluster heads only nodes
    S(i).type='N';
   
    temp_rnd0=i;
    %Random Election of Normal Nodes
    if (temp_rnd0>=m*n+1) 
        S(i).E=Eo;
        S(i).ENERGY=0;
        %%%%plot(S(i).xd,S(i).yd,'o');
        hold on;
    end
    %Random Election of Advanced Nodes
    if (temp_rnd0<m*n+1)  
        S(i).E=Eo*(1+a);
        E1=S(i).E;
        S(i).ENERGY=1;
        %%%%plot(S(i).xd,S(i).yd,'+');
         hold on;
    end
end

S(n+1).xd=sink.x;
S(n+1).yd=sink.y;
%%%%plot(S(n+1).xd,S(n+1).yd,'x');
    
        
%First Iteration
%estefigure(1);

%counter for CHs
countCHs=0;
%counter for CHs per round
rcountCHs=0;
cluster=1;

countCHs;
rcountCHs=rcountCHs+countCHs;
flag_first_dead=0;

for r=0:1:rmax
   % r

 %Election Probability for Normal Nodes
  pnrm=(p/(1+a*m*S(i).E));
  %Election Probability for Advanced Nodes
  padv= (p*(1+a)/(1+a*m*S(i).E));
    
  %Operation for heterogeneous epoch
  if(mod(r, round(1/pnrm) )==0)
    for i=1:1:n
        S(i).G=0;
        S(i).cl=0;
    end
  end

 %Operations for sub-epochs
 if(mod(r, round(1/padv) )==0)
    for i=1:1:n
        if(S(i).ENERGY==1)
            S(i).G=0;
            S(i).cl=0;
        end
    end
  end

 
hold off;

%Number of dead nodes
dead=0;
%Number of dead Advanced Nodes
dead_a=0;
%Number of dead Normal Nodes
dead_n=0;

%counter for bit transmitted to Bases Station and to Cluster Heads
packets_TO_BS=0;
packets_TO_CH=0;
%counter for bit transmitted to Bases Station and to Cluster Heads 
%per round
PACKETS_TO_CH(r+1)=0;
PACKETS_TO_BS(r+1)=0;

%estefigure(1);

for i=1:1:n
    %checking if there is a dead node
    if (S(i).E<=0)
 %este       plot(S(i).xd,S(i).yd,'red .');
        dead=dead+1;
        if(S(i).ENERGY==1)
            dead_a=dead_a+1;
        end
        if(S(i).ENERGY==0)
            dead_n=dead_n+1;
        end
        hold on;    
    end
    if S(i).E>0
        S(i).type='N';
        if (S(i).ENERGY==0)  
 %este       plot(S(i).xd,S(i).yd,'o');
        end
        if (S(i).ENERGY==1)  
  %este      plot(S(i).xd,S(i).yd,'+');
        end
        hold on;
    end
end
%este plot(S(n+1).xd,S(n+1).yd,'x');


STATISTICS(r+1).DEAD=dead;
DEAD(r+1)=dead;
DEAD_N(r+1)=dead_n;
DEAD_A(r+1)=dead_a;

%When the first node dies
if (dead==1)
    if(flag_first_dead==0)
        first_dead=r;
        flag_first_dead=1;
    end
end

countCHs=0;
cluster=1;
sum=0;
for i=1:1:n
if(S(i).E>0)
    sum=sum+S(i).E;
end
end
Eavg=sum/n;
for i=1:1:n
   if(S(i).E>0)
   temp_rand=rand;     
   if ( (S(i).G)<=0)

 %Election of Cluster Heads for normal nodes
 if( ( S(i).ENERGY==0 && ( temp_rand <= (( 3*pnrm / ( 1 - pnrm * mod(r,round(1/(3*pnrm))) ))*(Eavg/Eo)) ) )  )
     
     countCHs=countCHs+1;
     packets_TO_BS=packets_TO_BS+1;
     PACKETS_TO_BS(r+1)=packets_TO_BS;
     
     S(i).type='C';
     S(i).G=100;
     C(cluster).xd=S(i).xd;
     C(cluster).yd=S(i).yd;
 %este    plot(S(i).xd,S(i).yd,'k*');
     
     distance=sqrt( (S(i).xd-(S(n+1).xd) )^2 + (S(i).yd-(S(n+1).yd) )^2 );
     C(cluster).distance=distance;
     C(cluster).id=i;
     X(cluster)=S(i).xd;
     Y(cluster)=S(i).yd;
     cluster=cluster+1;
     
     %Calculation of Energy dissipated
     distance;
     if (distance>do)
         S(i).E=S(i).E- ( (ETX+EDA)*(4000) + Emp*4000*( distance*distance*distance*distance ));
     end
     if (distance<=do)
         S(i).E=S(i).E- ( (ETX+EDA)*(4000)  + Efs*4000*( distance * distance ));
     end
 end
 


 %Election of Cluster Heads for Advanced nodes
 if( ( S(i).ENERGY==1 && ( temp_rand <= (( 3*padv / ( 1 - padv * mod(r,round(1/(3*padv))) ))*(Eavg/E1)) ) )  )
        
            countCHs=countCHs+1;
            packets_TO_BS=packets_TO_BS+1;
            PACKETS_TO_BS(r+1)=packets_TO_BS;
            
            S(i).type='C';
            S(i).G=100;
            C(cluster).xd=S(i).xd;
            C(cluster).yd=S(i).yd;
 %este           plot(S(i).xd,S(i).yd,'k*');
            
            distance=sqrt( (S(i).xd-(S(n+1).xd) )^2 + (S(i).yd-(S(n+1).yd) )^2 );
            C(cluster).distance=distance;
            C(cluster).id=i;
            X(cluster)=S(i).xd;
            Y(cluster)=S(i).yd;
            cluster=cluster+1;
            
            %Calculation of Energy dissipated
            distance;
            if (distance>do)
                S(i).E=S(i).E- ( (ETX+EDA)*(4000) + Emp*4000*( distance*distance*distance*distance )); 
            end
            if (distance<=do)
                S(i).E=S(i).E- ( (ETX+EDA)*(4000)  + Efs*4000*( distance * distance )); 
            end
        end     
    
    end
  end 
end



STATISTICS(r+1).CLUSTERHEADS=cluster-1;
CLUSTERHS(r+1)=cluster-1;

%Election of Associated Cluster Head for Normal Nodes
for i=1:1:n
    if ( S(i).type=='N' && S(i).E>0 )
        if(cluster-1>=1)
            min_dis=sqrt( (S(i).xd-S(n+1).xd)^2 + (S(i).yd-S(n+1).yd)^2 );
            min_dis_cluster=1;
            for c=1:1:cluster-1
                temp=min(min_dis,sqrt( (S(i).xd-C(c).xd)^2 + (S(i).yd-C(c).yd)^2 ) );
                if ( temp<min_dis )
                    min_dis=temp;
                    min_dis_cluster=c;
                end
            end
            
            %Energy dissipated by associated Cluster Head
            min_dis;
            if (min_dis>do)
                S(i).E=S(i).E- ( ETX*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis));
            end
            if (min_dis<=do)
                S(i).E=S(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis));
            end
            %Energy dissipated
            if(min_dis>0)
                S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E- ( (ERX + EDA)*4000 );
                PACKETS_TO_CH(r+1)=n-dead-cluster+1;
            end
            
            S(i).min_dis=min_dis;
            S(i).min_dis_cluster=min_dis_cluster;
            
        end
    end
end
hold on;
%     if DEAD(r+1)==n
%         rmax=r+1;
%         break;
%     end
countCHs;
rcountCHs=rcountCHs+countCHs;
  sum=0;
for i=1:1:n
if(S(i).E>0)
    sum=sum+S(i).E;
end
end
avg=sum/n;
STATISTICS(r+1).AVG=avg;
uu(r+1).AVG=(STATISTICS(r+1).AVG)*100;
sum;
x5(r+1)=uu(r+1).AVG;

% if m==0.1
%     x=u;
% end
% if m==0.2
%     x1=u;
% end
% if m==0.3
%     x2=u;
% end

%Code for Voronoi Cells
%Unfortynately if there is a small
%number of cells, Matlab's voronoi
%procedure has some problems

%[vx,vy]=voronoi(X,Y);
%plot(X,Y,'r*',vx,vy,'b-');
% hold on;
% voronoi(X,Y);
% axis([0 xm 0 ym]);
   

end
alive4=u;
for i=1:1:rmax
    alive4(i)=n-DEAD(i);
%     if m==0.1
%         x(i)= alive(i);
%     end
% if m==0.2
%     x1(i)= alive(i);
% end
% if m==0.3
%     x2(i)= alive(i);
% end 
   
end
pack4=PACKETS_TO_BS;
save  pack4.mat;
save alive4.mat;
save  x5.mat;
% figure(3)
% r=1:1:rmax;
% plot(r,alive,'--');
% legend(['m=0.1','  ','\alpha=1']);
% ylabel('Number of AliveNodes')
% xlabel('Round Number')
% hold on;


%%

%%LEACH EEHC

clear all;
clear;

clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Field Dimensions - x and y maximum (in meters)
xm=100;
ym=100;

%x and y Coordinates of the Sink
sink.x=0.5*xm;
sink.y=0.5*ym;

%Number of Nodes in the field
n=100;

%Optimal Election Probability of a node
%to become cluster head
p=0;
%Energy Model (all values in Joules)
%Initial Energy 
Eo=0.1;
%Eelec=Etx=Erx
ETX=50*0.000000001;
ERX=50*0.000000001;
%Transmit Amplifier types
Efs=10*0.000000000001;
Emp=0.0013*0.000000000001;
%Data Aggregation Energy
EDA=5*0.000000001;

%Values for Hetereogeneity
%Percentage of nodes than are advanced
m=0.1;
%Percentage of nodes than are intermediate
x=0.2;
%\alpha
a=1;
%Beta
b=0.5;
%maximum number of rounds
rmax=500;
%variables
Ave_CH=0;
    sum=0;
count_ch=0;
Throughput=0;
 % packets size
Packet=4000;
%%%%%%%%%%%%%%%%%%%%%%%%% END OF PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%

%Computation of do
do=sqrt(Efs/Emp); 

%Creation of the random Sensor Network
figure(1);
 rand('seed',19)
for i=1:1:n
    S(i).xd=rand(1,1)*xm;
    XR(i)=S(i).xd;
    S(i).yd=rand(1,1)*ym;
    YR(i)=S(i).yd;
    S(i).G=0;
    S(i).E=0;
    %initially there are no cluster heads only nodes
    S(i).type='N';
   keep(i)=i;
    temp_rnd0=i;
    %Random Election of Normal Nodes
    if (temp_rnd0>=(x+m)*n+1) 
        S(i).E=Eo;
        S(i).ENERGY=0;
        %%%%plot(S(i).xd,S(i).yd,'o');
        hold on;
    end
      %Random Election of intermediate Nodes
    if (temp_rnd0<(m+x)*n+1) && (temp_rnd0>m*n)  
        S(i).E=Eo*(1+b);
        S(i).ENERGY=0.5;
        %%%%plot(S(i).xd,S(i).yd,'*');
         hold on;
    end
    %Random Election of Advanced Nodes
    if (temp_rnd0<m*n+1)  
        S(i).E=Eo*(1+a);
        S(i).ENERGY=1;
        %%%%plot(S(i).xd,S(i).yd,'D');
         hold on;
    end
   
end

S(n+1).xd=sink.x;
S(n+1).yd=sink.y;
%%plot(S(n+1).xd,S(n+1).yd,'x');
    
        
%First Iteration
figure(1);

%counter for CHs
countCHs=0;
%counter for CHs per round
rcountCHs=0;
cluster=1;
u=0;
countCHs;
rcountCHs=rcountCHs+countCHs;
flag_first_dead=0;
flag_first_Hdead=0;
flag_last_dead=0;
c=1;
  
for r=0:1:rmax
 distance=0.765*(xm/2);
K_opt=sqrt(n/(2*pi))*do*(xm/distance);
 p=K_opt/n;  
    
   for i=1:1:n
        if(S(i).E>0)
            holder(i)=S(i).E;
            id(i)=keep(i);
            node= struct('energy', holder, 'id',id);
            [energy,index] = sort([node.energy],'descend');  % Sort all energy values, largest first
        end
        
    end
    total=0;
    for k=1:length(node.energy)
        energy_level=sort(node.energy, 'descend');
        total=total + node.energy(k);
        
    end
       
        average=total/length(node.energy);
        
TEnergy(r+1)=total; 
AveEnergy(r+1)=average;
 
    r
  %Election Probability for Normal Nodes
  pnrm=( p/ (1+a*m+b*x) );
  %Election Probability for intermediate Nodes
  pint=( p*(1+b)/ (1+a*m+b*x) );
  %Election Probability for Advanced Nodes
  padv= ( p*(1+a)/(1+a*m+b*x) );
    
  %Operation for heterogeneous epoch
  if(mod(r, round(1/pnrm) )==0)
    for i=1:1:n
        S(i).G=0;
        S(i).cl=0;
    end
  end

 %Operations for sub-epochs
 if(mod(r, round(1/padv) )==0)
    for i=1:1:n
        if(S(i).ENERGY==1)
            S(i).G=0;
            S(i).cl=0;
        end
    end
 end

   %Operations for sub-epochs
 if(mod(r, round(1/pint) )==0)
    for i=1:1:n
        if(S(i).ENERGY==0.5)
            S(i).G=0;
            S(i).cl=0;
        end
    end
  end
 
hold off;
%Number of Half dead nodes
Hdead=0;
%Number of half dead Advanced Nodes
Hdead_a=0;
%Number of  half dead Normal Nodes
Hdead_n=0;
%Number of  half dead intermediate Nodes
Hdead_in=0;
%Number of dead nodes
dead=0;
%Number of dead Advanced Nodes
dead_a=0;
%Number of dead Normal Nodes
dead_n=0;
%Number of dead intermediate Nodes
dead_in=0;

%counter for bit transmitted to Bases Station and to Cluster Heads
packets_TO_BS=0;
packets_TO_CH=0;
%counter for bit transmitted to Bases Station and to Cluster Heads 
%per round
PACKETS_TO_CH(r+1)=0;
PACKETS_TO_BS(r+1)=0;

figure(1);

for i=1:1:n
    %Checking if the energy is less or reduced by half
    if (S(i).E<=(Eo/2)) && (S(i).E>0)
        plot(S(i).xd,S(i).yd,'yellow .');
        Hdead=Hdead+1;
        if(S(i).ENERGY==1)
            Hdead_a=Hdead_a+1;
        end
        if(S(i).ENERGY==0.5)
            Hdead_in=Hdead_in+1;
        end
        if(S(i).ENERGY==0)
            Hdead_n=Hdead_n+1;
        end
        hold on; 
    end   
    
    %checking the rate of energy dissipation in normal, intermediate and
    %advance nodes
   
   if (S(i).E<=Eo)||(S(i).E>Eo)
      if(S(i).ENERGY==0)
          RnEnergy(r+1)=S(i).E;
      end
      if (S(i).ENERGY==0.5)
          RINEnergy(r+1)=S(i).E;
     end
      if (S(i).ENERGY==1)
          RAEnergy(r+1)=S(i).E;
     end
    end

    %checking if there is a dead node
    if (S(i).E<=0)
        plot(S(i).xd,S(i).yd,'red .');
        dead=dead+1;
        if(S(i).ENERGY==1)
            dead_a=dead_a+1;
        end
        if(S(i).ENERGY==0)
            dead_n=dead_n+1;
        end
        if(S(i).ENERGY==0.5)
            dead_in=dead_in+1;
        end
        hold on;    
    end
    if S(i).E>0
        S(i).type='N';
        if (S(i).ENERGY==0)  
        plot(S(i).xd,S(i).yd,'o');
        end
        if (S(i).ENERGY==0.5)  
        plot(S(i).xd,S(i).yd,'p');
        end
        if (S(i).ENERGY==1)  
        plot(S(i).xd,S(i).yd,'D');
        end
        hold on;
    end
end
plot(S(n+1).xd,S(n+1).yd,'x');


HSTATISTICS(r+1).DEAD=Hdead;
HDEAD(r+1)=Hdead;
HDEAD_N(r+1)=Hdead_n;
HDEAD_IN(r+1)=Hdead_in;
HDEAD_A(r+1)=Hdead_a;

%When the first node is half dead
if (Hdead==1)
    if(flag_first_Hdead==0)
        first_Hdead=r;
        flag_first_Hdead=1;
    end
end

STATISTICS(r+1).DEAD=dead;
DEAD(r+1)=dead;
DEAD_N(r+1)=dead_n;
DEAD_IN(r+1)=dead_in;
DEAD_A(r+1)=dead_a;

%When the first node dies
if (dead==1)
    if(flag_first_dead==0)
        first_dead=r;
        flag_first_dead=1;
    end
end

%Number of  alive Nodes
alive=0;
%Number of  alive Normal Nodes
alive_n=0;
%Number of  alive intermediate Nodes
alive_in=0;
%Number of  alive advance Nodes
alive_a=0;



for i=1:1:n
    %checking number of alive node per round
    if (S(i).E>0)
        alive=alive+1;
        if(S(i).ENERGY==1)
            alive_a=alive_a+1;
        end
        if(S(i).ENERGY==0.5)
            alive_in=alive_in+1;
        end
        if(S(i).ENERGY==0)
            alive_n=alive_n+1;
        end
        hold on;    
    end
    %checking nodes status
if (S(i).E>0)
    nodes_status=1;
end
if (S(i).E<0)
    nodes_status=0;
end  
 STATISTICS(i).Status=nodes_status; 
 Status(i)=nodes_status;
        
    ASTATISTICS(r+1).Live=alive;
Live(r+1)=alive;
Live_n(r+1)=alive_n;
Live_in(r+1)=alive_in;
Live_a(r+1)=alive_a;
end
for i=1:1:n
%checking for last dead or alive node
if (alive==1 && S(i).E>0)
    if (S(i).ENERGY==1||S(i).ENERGY==0||S(i).ENERGY==0.5)
        plot(S(i).xd,S(i).yd,'green .');
        last_dead=r;
        Instability=last_dead-first_dead;
        flag_last_dead=1;
    end
    
end
end


countCHs=0;
cluster=1;
for i=1:1:n
   if(S(i).E>0)
   temp_rand=rand;     
   if ( (S(i).G)<=0)

if( ( S(i).ENERGY==0 && ( temp_rand <= ( pnrm / ( 1 - pnrm * mod(r,round(1/pnrm)) )) ) )  )

            countCHs=countCHs+1;
            packets_TO_BS=packets_TO_BS+1;
            PACKETS_TO_BS(r+1)=packets_TO_BS;
            
            S(i).type='C';
            S(i).G=100;
            C(cluster).xd=S(i).xd;
            C(cluster).yd=S(i).yd;
            plot(S(i).xd,S(i).yd,'k*');
            
%             distance=sqrt( (S(i).xd-(S(n+1).xd) )^2 + (S(i).yd-(S(n+1).yd) )^2 );
            C(cluster).distance=distance;
            C(cluster).id=i;
            X(cluster)=S(i).xd;
            Y(cluster)=S(i).yd;
            cluster=cluster+1;
            
            %Calculation of Energy dissipated
            distance;            
            if (distance>do)                                                
                S(i).E=S(i).E- ( ((n/cluster)-1)*(ETX+EDA)*(4000) +(n/cluster)*(EDA)*(4000)+(ETX+EDA)*(4000)  + Efs*4000*( distance^2));
            end
            if (distance<=do)
                S(i).E=S(i).E- ( ((n/cluster)-1)*(ETX+EDA)*(4000) +(n/cluster)*(EDA)*(4000)+(ETX+EDA)*(4000)  + Efs*4000*( distance)); 
            end
        end     
%     modular(r+1)=( 1 - 0.1 * mod(r,round(1/0.1)));
%     modular_ver(r+1)=( 1 - 0.1 * r);

 %Election of Cluster Heads for intermediate nodes
 if( ( S(i).ENERGY==0.5 && ( temp_rand <= ( pint / ( 1 - pint * mod(r,round(1/pint)) )) ) )  )
     
     countCHs=countCHs+1;
     packets_TO_BS=packets_TO_BS+1;
     PACKETS_TO_BS(r+1)=packets_TO_BS;
     
     S(i).type='C';
     S(i).G=100;
     C(cluster).xd=S(i).xd;
     C(cluster).yd=S(i).yd;
     plot(S(i).xd,S(i).yd,'k*');
     
%      distance=sqrt( (S(i).xd-(S(n+1).xd) )^2 + (S(i).yd-(S(n+1).yd) )^2 );
     C(cluster).distance=distance;
     C(cluster).id=i;
     X(cluster)=S(i).xd;
     Y(cluster)=S(i).yd;
     cluster=cluster+1;
     
     %Calculation of Energy dissipated
    distance;            
            if (distance>do)                                                
                S(i).E=S(i).E- ( ((n/cluster)-1)*(ETX+EDA)*(4000) +(n/cluster)*(EDA)*(4000)+(ETX+EDA)*(4000)  + Efs*4000*( distance^2));
            end
            if (distance<=do)
                S(i).E=S(i).E- ( ((n/cluster)-1)*(ETX+EDA)*(4000) +(n/cluster)*(EDA)*(4000)+(ETX+EDA)*(4000)  + Efs*4000*( distance)); 
            end
 end
    

 %Election of Cluster Heads for Advanced nodes
 if( ( S(i).ENERGY==1 && ( temp_rand <= ( padv / ( 1 - padv * mod(r,round(1/padv)) )) ) )  )
     
     countCHs=countCHs+1;
     packets_TO_BS=packets_TO_BS+1;
     PACKETS_TO_BS(r+1)=packets_TO_BS;
     
     S(i).type='C';
     S(i).G=100;
     C(cluster).xd=S(i).xd;
     C(cluster).yd=S(i).yd;
     plot(S(i).xd,S(i).yd,'k*');
     
%      distance=sqrt( (S(i).xd-(S(n+1).xd) )^2 + (S(i).yd-(S(n+1).yd) )^2 );
     C(cluster).distance=distance;
     C(cluster).id=i;
     X(cluster)=S(i).xd;
     Y(cluster)=S(i).yd;
     cluster=cluster+1;
     
     %Calculation of Energy dissipated
     distance;            
            if (distance>do)                                                
                S(i).E=S(i).E- ( ((n/cluster)-1)*(ETX+EDA)*(4000) +(n/cluster)*(EDA)*(4000)+(ETX+EDA)*(4000)  + Efs*4000*( distance^2));
            end
            if (distance<=do)
                S(i).E=S(i).E- ( ((n/cluster)-1)*(ETX+EDA)*(4000) +(n/cluster)*(EDA)*(4000)+(ETX+EDA)*(4000)  + Efs*4000*( distance)); 
            end
 end
    
    end
  end 
end

%Checking average number of ClusterHeads per epoch    

sum=sum+cluster;
count_ch=count_ch+1;

if  count_ch==10
    Ave_CH=(sum*0.1)/(1+(m*a)+(b*x));
    Throughput=Ave_CH*4;
    
    STATISTICS(r+1).ave_clustHd=Ave_CH;
    ave_ch(r+1)=Ave_CH;
    STATISTICS(r+1).throughput=Throughput;
    Clust_throughput(r+1)=Throughput;
    Ave_CH=0;
    sum=0;
    count_ch=0;
    
end

STATISTICS(r+1).CLUSTERHEADS=cluster-1;
CLUSTERHS(r+1)=cluster-1;
countmember=0;

%Election of Associated Cluster Head for Normal Nodes
for i=1:1:n
    if ( S(i).type=='N' && S(i).E>0 )
        if(cluster-1>=1)
            min_dis=sqrt( (S(i).xd-S(n+1).xd)^2 + (S(i).yd-S(n+1).yd)^2 );
            min_dis_cluster=1;
            for c=1:1:cluster-1
                temp=min(min_dis,sqrt( (S(i).xd-C(c).xd)^2 + (S(i).yd-C(c).yd)^2 ) );
                if ( temp<min_dis )
                    min_dis=temp;
                    min_dis_cluster=c;
                end
            end
            
            %Energy dissipated by associated Cluster Head
            min_dis;
            min_dis=(xm*ym)/(2*(cluster-1)*pi);
            if (min_dis>do)
                S(i).E=S(i).E- ( ETX*(4000) + Emp*4000*( min_dis * min_dis ));
            end
            if (min_dis<=do)
                S(i).E=S(i).E- ( ETX*(4000) + Efs*4000*( min_dis));
            end
            %Energy dissipated by associated Cluster Head
            min_dis;
            if (min_dis>do)
                S(i).E=S(i).E- ( ETX*(Packet) + Emp*Packet*( min_dis * min_dis * min_dis * min_dis));
                Energy_member=ETX*(Packet) + Emp*Packet*( min_dis * min_dis * min_dis * min_dis);
            end
            if (min_dis<=do)
                S(i).E=S(i).E- ( ETX*(Packet) + Efs*Packet*( min_dis * min_dis));
                Energy_member=ETX*(Packet) + Emp*Packet*( min_dis * min_dis);
            end
            %Energy dissipated
            if(min_dis>0)
                S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E- ( (ERX + EDA)*Packet );
                PACKETS_TO_CH(r+1)=n-dead-cluster+1;
            end
            
            S(i).min_dis=min_dis;
            S(i).min_dis_cluster=min_dis_cluster;
            
        end
    end
end
hold on;

countCHs;
rcountCHs=rcountCHs+countCHs;

  sum=0;
for i=1:1:n
if(S(i).E>0)
    sum=sum+S(i).E;
end
end
avg=sum/n;
STATISTICS(r+1).AVG=avg;
uu(r+1).AVG=(STATISTICS(r+1).AVG)*100;
sum;
alive5=u;
x6(r+1)=uu(r+1).AVG;

% if m==0.1
%     x=u;
% end
% if m==0.2
%     x1=u;
% end
% if m==0.3
%     x2=u;
% end

%Code for Voronoi Cells
%Unfortynately if there is a small
%number of cells, Matlab's voronoi
%procedure has some problems

%[vx,vy]=voronoi(X,Y);
%plot(X,Y,'r*',vx,vy,'b-');
% hold on;
% voronoi(X,Y);
% axis([0 xm 0 ym]);
   

end
for i=1:1:rmax
    alive5(i)=n-DEAD(i);
%     if m==0.1
%         x(i)= alive(i);
%     end
% if m==0.2
%     x1(i)= alive(i);
% end
% if m==0.3
%     x2(i)= alive(i);
% end 
   
end
pack5=PACKETS_TO_BS;
save  pack5.mat;
save alive5.mat;
save  x6.mat;
% figure(3)
% r=1:1:rmax;
% plot(r,alive,'--');
% legend(['m=0.1','  ','\alpha=1']);
% ylabel('Number of AliveNodes')
% xlabel('Round Number')
% hold on;


%%
clear all;
clc;
load ('alive.mat');
load ('alive1.mat');
load ('alive2.mat');
%load ('alive3.mat');
load ('alive4.mat');
load ('alive5.mat');
load ('x1.mat');
load ('x2.mat');
load ('x3.mat');
%load ('x4.mat');
load ('x5.mat');
load ('x6.mat');
load ('pack.mat');
load ('pack2.mat');
load ('pack1.mat');
%load ('pack3.mat');
load ('pack4.mat');
load ('pack5.mat');

figure(3)
r=1:1:rmax;
plot(r,alive,'--',r,alive1,'g',r,alive2,':',r,alive4,'--+',r,alive5,'--*');
legend('P-SEP','SEP','LEACH-DCHS','M-SEP','LEACH-EEHC');
%,'A-LEACH' %,r,alive3,'--o'
ylabel('Number of AliveNodes');
xlabel('Round Number');
hold on;

figure(4)
r=1:1:rmax;
plot(x1(r),r,'--',x2(r),r,'g',x3(r),r,':',x5(r),r,'--+',x6(r),r,'--*');
legend('P-SEP','SEP','LEACH-DCHS','M-SEP','LEACH-EEHC');
% ,'A-LEACH'   %,x4(r),r,'--o'
xlabel('Average Energy of Each Node');
ylabel('Round Number');   
hold on;

figure(5)
r=1:1:rmax;
plot(r,pack(r),'--',r,pack1(r),'g',r,pack2(r),':',r,pack4(r),'--+',r,pack5(r),'--*');
legend('P-SEP','SEP','LEACH-DCHS','M-SEP','LEACH-EEHC');  
% ,'A-LEACH'  %,r,pack3(r),'--o'
ylabel('Packets to Sink');
xlabel('Round Number');
hold on;

%%

clear all;
clc;
load ('alive5.mat');


figure(3)
r=1:1:rmax;
plot(r,alive5,'--');
legend('SEP');
%,'A-LEACH' %,r,alive3,'--o'
ylabel('Number of AliveNodes');
xlabel('Round Number');
hold on;


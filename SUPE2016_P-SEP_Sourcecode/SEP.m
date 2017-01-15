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
 Eo           =       0.5;          %Initial Energy
ETX          =        50*(10^-9);   %Eelec=Etx=Erx
ERX          =        50*(10^-9);   %Eelec=Etx=Erx
Efs          =        10*(10^-12);  %Transmit Amplifier types
Emp          =        0.0013*(10^-12);
EDA          =        5*(10^-9);    %Data Aggregation Energy
a             =       1;            %\alpha - Percentage of nodes than are advanced - Values for Hetereogeneity
rmax         =        3000;         %maximum number of rounds
z            =        1;
do           =        sqrt(Efs/Emp);%distancia
rp           =        m*n+1;        %probabilidad de nodes advanced
j            =        0;            %variable para ordenar nodos advanced
u            =        0;
 x1=0;
% %Thresholod for transmiting data to the cluster head
% t_h          =        10;           %%%%%%Hard Thres%%%%hold H(t)
% t_s          =        5;             %%%%%%Soft thres%%%%hold  S(t)
% sv           =        0;             %%%%%%previously Sensed value S(v)
% %methane level range
% l_r1         =        2.5;
% l_r          =        12.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%  END OF PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%  Creation of the random Sensor Network  %%%%%%%%%%%%%%
figure(1);
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
        S(i).ENERGY=0;
        S(i).t=0;
        %%%%plot(S(i).xd,S(i).yd,'o');
        hold on;
    end
    
    %Random Election of Advanced Nodes
    if (temp_rnd0<m*n+1)
        S(i).E=Eo*(1+a)
        S(i).ENERGY=1;
        S(i).t=0;
        %%%%plot(S(i).xd,S(i).yd,'+');
        if (mod(i,j)>=1) && (mod(i,j)<l)
            S(i).xd=rpr1;
            S(i).yd=distancia1(z);
            %%%%plot(S(i).xd,S(i).yd,'+');
            z=z+1;
        end
        
        if (mod(i,j)==0)
            S(i).xd=rpr1;
            S(i).yd=distancia1(z);
            z=1;
            d=d+1;
            rpr1=rpr*d;
        end
    end
end



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
inisv=0;
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
    %Election Probability for Normal Nodes
    pnrm=( p/ (1+a*m) );
    %Election Probability for Advanced Nodes
    padv= ( p*(1+a)/(1+a*m) );
    
    %Operation for heterogeneous epoch
    
     if mod(r,2)==0
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
    packets_TO_BS1=0;
    packets_TO_CH=0;
    %counter for bit transmitted to Bases Station and to Cluster Heads
    %per round
    PACKETS_TO_CH(r+1)=0;
    PACKETS_TO_BS(r+1)=0;
    PACKETS_TO_BS1(r+1)=0;
    
    figure(1);
    
    for i=1:1:n
        %checking if there is a dead node
        if (S(i).E<=0)
            plot(S(i).xd,S(i).yd,'red .');
            dead=dead+1;
            if(S(i).ENERGY==1)   %%  si es =1 quiere decir q se mueve es un nodo adv
                dead_a=dead_a+1;
            end
            if(S(i).ENERGY==0) %%  si es =0 quiere decir q se no se mueve es un nodo normal
                dead_n=dead_n+1;
            end
            hold on;
        end
        if S(i).E>0
            S(i).type='N';   %%%% quiere decir q los nodos aun no estan muertos
            if (S(i).ENERGY==0)   %%  dibuja los nodos normales con energia 0(identificacion)
                plot(S(i).xd,S(i).yd,'o','MarkerEdgeColor','b',...
                    'MarkerFaceColor','w',...
                    'MarkerSize',5);
            end
            if (S(i).ENERGY==1)
                plot(S(i).xd,S(i).yd,'D','MarkerEdgeColor','g',...
                    'MarkerFaceColor','g',...
                    'MarkerSize',6);      
                %             circle(S(i).xd,S(i).yd,radio)
                hold on;
            end
        end
    end
    plot(S(n+1).xd,S(n+1).yd,'-k*','MarkerSize',15,'MarkerEdgeColor','r');  %%  dibuja el nodo  sink
    
    
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
    if(S(i).E>0.001)
        temp_rand=rand;
        if ( (S(i).G)==0)
            %Selecciona-Cluster Heads of normal nodes and advance nodes
            if( ( S(i).ENERGY==0 && ( temp_rand <= ( pnrm / ( 1 - pnrm * mod(r,round(1/pnrm)) )) ) )  ) || ( ( S(i).ENERGY==1 && ( temp_rand <= ( padv / ( 1 - padv * mod(r,round(1/padv)) )) ) )  )
                W(i).x=S(i).xd;
                W(i).y=S(i).yd;
%                 if (mod(r,2)==0)
for j=1:1:i
    distanciaz=sqrt( (S(i).xd-(W(j).x) )^2 + (S(i).yd-(W(j).y) )^2 );
    if ((distanciaz>5))%&((distanciaz<80))
        countCHs=countCHs+1;
        packets_TO_BS=packets_TO_BS+1;
        PACKETS_TO_BS(r+1)=packets_TO_BS;
        S(i).type='C';
        S(i).G=100;
        C(cluster).xd=S(i).xd;
        C(cluster).yd=S(i).yd;
        plot(S(i).xd,S(i).yd,'k*','MarkerEdgeColor','red',...
            'MarkerFaceColor','red',...
            'MarkerSize',6);
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
         break;
        
        %                         end
    end
end
            end
        end
    end
end


    
    STATISTICS(r+1).CLUSTERHEADS=cluster-1;
    CLUSTERHS(r+1)=cluster-1;
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
%                  if (sv==0)                   
                    if (min_dis>do)
                        S(i).E=S(i).E-( ETX*(4000) + Emp*4000*( min_dis^4));
                    end
                    if (min_dis<=do)
                        S(i).E=S(i).E-( ETX*(4000) + Efs*4000*( min_dis^2));
                    end
%                  end
              
                %Energy dissipated
%                  if (sv==0)
                    if(min_dis>0)
                        S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E- ((ERX + EDA)*4000 );
                        PACKETS_TO_CH(r+1)=n-dead-cluster+1;
                   end
                    
%                 end
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
end


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

save alive.mat;
save  x1.mat;

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

clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Field Dimensions - x and y maximum (in meters)
xm=100;
ym=100;

%x and y Coordinates of the Sink
sink.x=0.5*xm;
sink.y=0.5*ym;

%Number of Nodes in the field
n=100
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
rmax=3000;

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
        %%%%plot(S(i).xd,S(i).yd,'o');
        hold on;
    end
    %Random Election of Advanced Nodes
    if (temp_rnd0<m*n+1)  
        S(i).E=Eo*(1+a)
        S(i).ENERGY=1;
        %%%%plot(S(i).xd,S(i).yd,'+');
         hold on;
    end
end

S(n+1).xd=sink.x;
S(n+1).yd=sink.y;
%%%%plot(S(n+1).xd,S(n+1).yd,'x');
    
        
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

figure(1);

for i=1:1:n
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
        hold on;    
    end
    if S(i).E>0
        S(i).type='N';
        if (S(i).ENERGY==0)  
        plot(S(i).xd,S(i).yd,'o');
        end
        if (S(i).ENERGY==1)  
        plot(S(i).xd,S(i).yd,'+');
        end
        hold on;
    end
end
plot(S(n+1).xd,S(n+1).yd,'x');


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
     plot(S(i).xd,S(i).yd,'k*');
     
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
            plot(S(i).xd,S(i).yd,'k*');
            
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
%%LEACH

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
m=0.2;
%\alpha
a=1;
u=0;
x3=0;
%maximum number of rounds
rmax=3000;

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
        plot(S(i).xd,S(i).yd,'o');
        hold on;
    end
    %Random Election of Advanced Nodes
    if (temp_rnd0<m*n+1)  
        S(i).E=Eo*(1+a)
        S(i).ENERGY=1;
        plot(S(i).xd,S(i).yd,'+');
        hold on;
    end
end

S(n+1).xd=sink.x;
S(n+1).yd=sink.y;
plot(S(n+1).xd,S(n+1).yd,'x');
    
        
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
        plot(S(i).xd,S(i).yd,'red .');
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
        plot(S(i).xd,S(i).yd,'o');
        end
        if (S(i).ENERGY==1)  
        plot(S(i).xd,S(i).yd,'+');
        end
        hold on;
    end
end
plot(S(n+1).xd,S(n+1).yd,'x');


STATISTICS(r+1).DEAD=dead;
DEAD(r+1)=dead;
DEAD_N(r+1)=dead_n;
DEAD_A(r+1)=dead_a;

%When the first node dies
if (dead==1)
    if(flag_first_dead==0)
        first_dead=r
        flag_first_dead=1;
    end
end

countCHs=0;
cluster=1;
for i=1:1:n
   if(S(i).E>0)
   temp_rand=rand;     
   if ( (S(i).G)<=0)

 %Election of Cluster-Heads
 if(temp_rand<= (p/(1-p*mod(r,round(1/p)))))
            countCHs=countCHs+1;
            packets_TO_BS=packets_TO_BS+1;
            PACKETS_TO_BS(r+1)=packets_TO_BS;
            
            S(i).type='C';
            S(i).G=round(1/p)-1;
            C(cluster).xd=S(i).xd;
            C(cluster).yd=S(i).yd;
            plot(S(i).xd,S(i).yd,'k*');
            
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
 
save('alive2');
save('x3');
% figure(5)
% r= 1:1:rmax;
% plot(r,alive,'r');  
%%

load ('alive.mat');
load alive1.mat;
load alive2.mat;
load x1.mat;
load x2.mat;
load x3.mat;

figure(3)
r=1:rmax;
plot(r,alive,'--',r,alive1,'g',alive2,r,':'),
legend('P-SEP','SEP','LEACH');
ylabel('Number of AliveNodes')
xlabel('Round Number')
hold on;

figure(4)
r=1:1:rmax;
plot(x1(r),r,'--',x2(r),r,'g',x3(r),r,':');
legend('P-SEP','SEP','LEACH');
    xlabel('Average Energy of Each Node');
    ylabel('Round Number');   
     hold on;




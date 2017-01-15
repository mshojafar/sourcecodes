clear;
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  %Field Dimensions - x and y maximum (in meters)
xm           =        200;
ym           =        200;
axis([0 xm 0 ym]);
grid on;
for loop = 1:3
   if loop==1
       m=0.1;
   end
    if loop==2
       m=0.2;
     end
    if loop==3
       m=0.3;
   end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% VARIABLES  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n            =        100;          %Number of Nodes in the field
radio        =        n/10;         %Radio
% m            =        0.1;          %Probbilidad de que un nodo se convierta en Cluster Head
p            =        0.1;          %Probabilidad
Eo           =        0.5;          %Initial Energy
ETX          =        50*(10^-9);   %Eelec=Etx=Erx
ERX          =        50*(10^-9);   %Eelec=Etx=Erx
Efs          =        10*(10^-12);  %Transmit Amplifier types
Emp          =        0.0013*(10^-12);
EDA          =        5*(10^-9);    %Data Aggregation Energy
a            =        3;            %\alpha - Percentage of nodes than are advanced - Values for Hetereogeneity
rmax         =        7000;         %maximum number of rounds
z            =        1;
do           =        sqrt(Efs/Emp);%distancia
rp           =        m*n+1;        %probabilidad de nodes advanced
j            =        0;            %variable para ordenar nodos advanced
u            =        0;
ENERGIAT     =        50;
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
sink.x=0.5*xm;
sink.y=0.5*ym;
S(n+1).xd=sink.x;
S(n+1).yd=sink.y;

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
            if (S(i).ENERGY==0)   %%  dibuja los nodos normales ocn energia 0
                  S(i).xd=rand(1,1)*xm;
                  S(i).yd=rand(1,1)*ym;
                  plot(S(i).xd,S(i).yd,'o','MarkerEdgeColor','b',...
                    'MarkerFaceColor','w',...
                    'MarkerSize',5);
            end
            if (S(i).ENERGY==1)
                plot(S(i).xd,S(i).yd,'D','MarkerEdgeColor','g',...
                    'MarkerFaceColor','g',...
                    'MarkerSize',6);       %%  dibuja los nodos adv ocn energia 1
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
               %Election of Cluster Heads for all the nodes
                if( ( S(i).ENERGY==0 && ( temp_rand <= ( pnrm / ( 1 - pnrm * mod(r,round(1/pnrm)))))))|| ((S(i).ENERGY==1 && ( temp_rand <= ( padv / ( 1 - padv * mod(r,round(1/padv)))))))
                    
                    countCHs=countCHs+1;
                    packets_TO_BS=packets_TO_BS+1;
                    PACKETS_TO_BS(r+1)=packets_TO_BS;
                                                                             
                    S(i).type='C';
                    S(i).G=100;
                    C(cluster).xd=S(i).xd;
                    C(cluster).yd=S(i).yd;
                    plot(S(i).xd,S(i).yd,'k*','MarkerEdgeColor','black',...
                    'MarkerFaceColor','black',...
                    'MarkerSize',6);
                         
                
                
                    distance=sqrt( (S(i).xd-(S(n+1).xd) )^2 + (S(i).yd-(S(n+1).yd) )^2 );
                    C(cluster).distance=distance;
                    C(cluster).id=i;
                    X(cluster)=S(i).xd;
                    Y(cluster)=S(i).yd;
                    cluster=cluster+1;
                    
                    %Calculation of Energy dissipated
                    distance;
                    if (distance>do)
                        S(i).E=S(i).E-((ETX+EDA)*(4000)+ Emp*4000*( distance^4));
                    end
                    if (distance<=do)
                        S(i).E=S(i).E-((ETX+EDA)*(4000)+ Efs*4000*( distance^2 ));
                    end
                end                             
            end
        end
        energiat=S(i).E;
    end
    ENERGIAT(r+1)=energiat*100;
    
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
           energiat=S(i).E*100;
    end
 ENERGIAT(r+1)=energiat;
    if DEAD(r+1)==n
        rmax=r+1;
        break;
    end
    
    hold on;
    
    countCHs;
    rcountCHs=rcountCHs+countCHs;
    
    % axis([0 100 0 100]);
    grid on;
 end


alive=u;
energiar=u;
for i=1:1:rmax
    alive(i)=n-DEAD(i); 
   energiar(i)=ENERGIAT(i);
end
 


figure(2)
r=1:1:rmax;
% plot(energiar,r,'color',rand(1,3));
plot(r,alive,'color',rand(1,3));
hold on;
xlabel('Number of rounds')
ylabel('Number of Alive Nodes')
legend(['m=' num2str(m),' ','\alpha=' num2str(a)]);
grid on;
end

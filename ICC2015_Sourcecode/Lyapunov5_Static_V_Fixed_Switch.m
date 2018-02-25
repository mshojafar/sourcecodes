%variables inizialization utilizzo flessibile delle code , non ha vincoli
%di tempo di servizio.% CHE THROUPUT GARANTISCONO? E COSTO AL VARIARE v
clc
clear all
close all
cc=0;
N = 1;                                                                    %applications pag 479
M = 100;                                                                   %Servers pag 479
Ts = 2;                                                                    %time slot (sec)
TT = 2000;
T=2000;
k_e=0.005; %(milliJ/(MHz)^2)
%Lamda = 2000.*ones(N,TT);                                                  %Number of slot
%A = 2.*Lamda.*rand(N,TT);                                                  %Arrival job uniformed distribuited pag 481
%A=4000.*rand(N,TT);                                                        %PMR = 1
%A = 3000 +(4000-2000).*rand(N,TT);                                         %PMR = 1.25
%A = 2000 +(6000-2000).*rand(N,TT) ;                                        %PMR = 1.5
% A = 1000 +(7000-1000).*rand(N,TT);                                        %PMR = 1.75
%A = 8000.*rand(N,TT);                                                        %PMR = 2
%save('Avalue4000PMR2.mat','A');
%A=load('Avalue4000PMR2.mat','A');

%L_tot_vect = 6+ round((10-6) * rand(N,TT));                         % PMR=1.25, L_tot=8+-2
%save('workloads8PMR125No2000.mat', 'L_tot_vect');                  %loading 100 fixed random workloads
A=load('workloads8PMR125No2000.mat','L_tot_vect');
A=A.L_tot_vect;


%F1 = [1.6*10^9 :0.2*10^9:2.6*10^9];                                         %frequecny range for every VM allocated onto the j-server
F1 = [5,50,70,90,105];                                         %( Gbyte/s)frequecny range for every VM allocated onto the j-server
I = ones(M,TT);                                                            %resource decision allocation
P_cpu1 = zeros(N,M,length(F1));
mu1=zeros(N,M,length(F1)) ;                                                  %power consuption for each VM onto each j server
P_min1 = 0.*ones(N,M,length(F1));                                          %minimum power for the considered VM allocated onto the j-server
P_max1 = 60.*ones(N,M,length(F1));                                          %maximum power for the considered VM allocated onto the j-server

MU=zeros(N,M,TT);                                                          %service rate (request/slots)
ww = zeros(1,N);
AU = zeros(N,M);
AAU = zeros(1,N);
SAU = zeros(1,M);
AAS = zeros(1,N);
Avg_ADELAY = zeros(1,N);
avg_MU=zeros(N,M);
avg_AMU=zeros(1,N);
avg_SMU=zeros(1,M);
Avg_SDELAY = zeros(1,M);
Avg_SDelay=zeros(1,M);

KU=  zeros(N,M);
Delay=zeros(1,N);
Delay2=zeros(N,M);
%Avg_Delay1=0;
%UUTILITY = 0;

alpha = ones(1,N);                                                         %throughput utility weights
%V1 = [1:100:1000];                                                         %control parameter for DCA
V = 100;
%V=length(V1)
beta = 1;                                                                  % non-negative normalizing weight
W = zeros(N,TT);                                                           %Buffer Dimension
K = zeros(N,M,TT);                                                         % ==R(i,j) number of requests for application i that are routed from the R(i) router Buffer to the j-server in slot t.
U = zeros(N,M,TT);                                                          %queing dynamics for the request of application i at server j
P = zeros(M,TT);                                                            % Switch consumed by each server in each time slot for s
Switch=zeros(M,TT);                                                            % power consumed by each server in each time slot for s
r = zeros(1,N);                                                            %average expected rate of admitted request for apllication i
e = zeros(1,M);                                                            %the time average expected power consuption of server j
SW= zeros(1,M);                                                            %the time average expected power consuption of server j
R = zeros(N,TT);                                                           %the time average expected switch power consuption of server j
Freq1 = zeros(N,M,length(F1));
Freqtemp = zeros(M,TT);
temp1 = zeros(N,M,length(F1));
active_servers = zeros(N+1,M,T);
M1=0;
M2=0;
M3=0;
l=0;
on_servers_list=zeros(1,M);
on_servers = zeros(N+1,M);                                                 % on-server in each different frame T
off_servers = zeros(N+1,M);                                                % off-server in each different frame T
sleeping_servers= zeros(N+1,M);                                            % Hybernated-server in each different frame T
 a = rand(1,M);                                                             % applications variable indicator
 save('avalue.mat','a');
a=load('avalue.mat','a');
a=a.a;
a(:)=(a(:) >= 0.4);
ACTIVE_Num=1;
a1=zeros(1, M);
%a1(:) = a(:);
count=0;                                                                   % static servers 1: static server is and 0: this servers is Dynamic
for j=1:M
    if (a(j)==1) && (count<ACTIVE_Num)
        a1(j)=1;
        count=count+1;
    end
end
for j=1:M
    if (a(j)==1)
        on_servers(1,j)=j;
        l=l+1;
        on_servers_list(j)=j;
    end
end

a_temp = zeros(1,M);                                                        %  temp applications variable indicator
%a_temp(:)=(a_temp(:) >= 0.6);

for i = 1:N
    for j = 1:M
        Freq1(i,j,:)= F1(:);
    end
end

for i=1:N
    for j=1:M
        mu1(i,j,:) = (8*10^(-3).*F1(:)) +0.76;                  % servers queue for each application service rate (5,0.8)-----(105,1.6)
        P_cpu1(i,j,:)=P_min1(i,j,:) + (50*10^(-4)).*((Freq1(i,j,:)-5).^2);     % CPU power of each virtual machine among each j-server pag 480
    end
end

i=0;
j=0;
F3= [0, 5,50,70,90,105];    %frequecny range for every VM allocated onto the j-server
temp3 = zeros(N,M,length(F3));
P_cpu3 = zeros(N,M,length(F3));                                            % power consuption for each VM onto each j server
P_min3 = 0.*ones(N,M,length(F3));                                          %minimum power for the considered VM allocated onto the j-server
P_max3 = 60.*ones(N,M,length(F3));                                          %maximum power for the considered VM allocated onto the j-server
mu3 = zeros(N,M,length(F3))  ;
Freq3 = zeros(N,M,length(F3));
for i = 1:N
    for j = 1:M
        Freq3(i,j,:)= F3(:);
    end
end

for i=1:N
    for j=1:M
        mu3(i,j,2:1:length(F3)) = (8*10^(-3).*F3(2:1:length(F3)))+0.76; % servers queue for each application service rate (5,0.8)-----(105,1.6)
        mu3(i,j,1) = 0;
        P_cpu3(i,j,2:1:length(F3))=P_min3(i,j,2:1:length(F3)) + (50*10^(-4)).*((Freq3(i,j,2:1:length(F3))-(5)).^2);  % CPU power of each virtual machine among each j-server pag 480
        P_cpu3(i,j,1)=0;
        
    end
end
F2= [0, 105];    %frequecny range for every VM allocated onto the j-server (STATIC server)
P_cpu2 = zeros(N,M,length(F2));
P_min2 = 0.*ones(N,M,length(F2));                                          %minimum power for the considered VM allocated onto the j-server
P_max2 = 60.*ones(N,M,length(F2));                                          %maximum power for the considered VM allocated onto the j-server
mu2 = zeros(N,M,length(F2))  ;
%temp2 = zeros(N,M,length(F2));
Freq2 = zeros(N,M,length(F2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%(STATIC server)
for i=1:N
    for j=1:M
        mu2(i,j,2) = (8*10^(-3).*F2(2))+0.76;                                         % servers queue for each application service rate (5,0.8)-----(105,1.6)
        mu2(i,j,1) = 0;
        P_cpu2(i,j,2)=P_min2(i,j,length(F2)) + (50*10^(-4)).*((Freq2(i,j,length(F2))-(5)).^2);  % CPU power of each virtual machine among each j-server pag 480
        P_cpu2(i,j,1)=0;
        
    end
end

i=0;
j=0;
UU=zeros(N,M,TT);
UU1=zeros(N,M,TT);
B = [1:1:TT/T];
L_off=zeros(1,length(B));
MM1=zeros(1,length(B));
MM2=zeros(1,length(B));
MM3=zeros(1,length(B));
bb=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%main body code
for t = 1 : TT
    if (mod(t,T)~= 0)
        UU=zeros(N,M,TT);
        for jtest=1:M
            if (a(jtest)==1)
                UU(:,jtest,t)=U(:,jtest,t);
            else
                UU(:,jtest,t)=nan;
            end
        end
        %if (t>1)
        %UU(~UU)=nan;
        %end
        for i= 1 : N
            
            %
            %j
            t
            if(W(i,t) > V.*alpha(i))
                R(i,t) = 0;
                
            else
                R(i,t) = A(i,t);
            end
            
            %list U of active (ON) servers.
            
            
            [U1, indexu1]=min(UU(i,:,t));% page 483 routing
            if (W(i,t) > U1) || (isnan(U1))
                K(i,indexu1,t) = W(i,t);
            else
                K(i,:,t)=0;
            end
            
            UU(:,indexu1,t)=nan;
            
            %                 [U1, ind]=min(U(i,:,t));% page 483 routing
            %                 if (W(i,t) > U1)
            %                     for ikl=1:M
            %                         for mmmm=1:N
            %                             K(mmmm,ikl,t)=K(mmmm,ikl,t-1);
            %                         end
            %                     end
            %                     K(i,ind,t) = W(i,t);
            %
            %                 else
            %                     K(i,:,t)=0;
            %                 end
            
            
            %a(i,j) = (a(i,j) > 0.5); % garantee to set more than half to be -on
            % fix number of active server, they are not M , they are
            % obviously less
            
            %if sum(a(i,j)*K(i,j,t)<=W(i,t))
            %if ((W(i,t)-sum(K(i,:,t)))>=0)
            
            % W(i,t+1) = W(i,t) - sum(a(i,j).*K(i,:,t)) + R(i,t);
            W(i,t+1) = W(i,t) - (a*K(i,:,t)') + R(i,t);                %% Backlog Queue 1    (2) pag 481 consider only active server
            % else
            %    W(i,t+1)=R(i,t);
            % end
            
        end % for N-application
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% after routing inside each server
        for j= 1 : M
            if (a(j)==1) && (a1(j)==0)                                      % for Dynamic ON server
                [U_max1 , index_U_max1]=max(U(:,j,t));                   % 0=<index_U_max<=N
                
                for zzz=1:length(F1)
                    temp1(index_U_max1,j,zzz) = U_max1.*(min(mu1(index_U_max1,j,zzz),U_max1)) - V.*beta.*(P_cpu1(index_U_max1,j,zzz));                  % Resource allocation pag 483
                    
                end
                
                [value1 , index1]=max(temp1(index_U_max1,j,:));            % 1=<index<=length(F)
                
                MU(index_U_max1,j,t) = min(mu1(index_U_max1,j,index1),U_max1);
                P(j,t)= P_cpu1(index_U_max1,j,index1);
                if (t==1)
                    Switch(j,t)=k_e*(F1(index1))^2;
                    Freqtemp(j,2)=F1(index1);
                else
                    Switch(j,t)=k_e*((F1(index1)-Freqtemp(j,t-1)))^2;
                    Freqtemp(j,t)=F1(index1);
                end
                for i= 1 : N
                    if (i==index_U_max1)
                        U(index_U_max1,j,t+1) = max((U(index_U_max1,j,t) - MU(index_U_max1,j,t)),0) + K(index_U_max1,j,t);      %% Backlog Queue 2    (5) pag 481
                    else
                        U(i,j,t+1) = max(U(i,j,t),0) + K(i,j,t);           % Backlog Queue 2    (5) pag 481
                    end
                end
                % U(index_U_max,j,t+1) = max(U(index_U_max,j,t) - MU(index_U_max,j,t),0) + a(i,j).*K(index_U_max,j,t);      %% Backlog Queue 2    (5) pag 481
                
   
                
       % for STATIC server % for STATIC server % for STATIC server % for STATIC server % for STATIC server % for STATIC server
              
                
            elseif  (a(j)==1) && (a1(j)==1)                                  % for Static ON server
                
                [U_max2 , index_U_max2]=max(U(:,j,t));                  % 0=<index_U_max<=N
                
                MU(index_U_max2,j,t) = min(mu2(index_U_max2,j,length(F2)),U_max2);
                %P(j,t)= P_cpu2(index_U_max2,j,length(F2));
                P(j,t)= P_max2(index_U_max2,j,length(F2));
                if (t==1)
                    Switch(j,t)=k_e*(F2(length(F2)))^2;
                    Freqtemp(j,2)=F2(length(F2));
                else
                    Switch(j,t)=k_e*((F2(length(F2))-Freqtemp(j,t-1)))^2;
                    Freqtemp(j,t)=F2(length(F2));
                end
                for i= 1 : N
                    if (i==index_U_max2)
                        U(index_U_max2,j,t+1) = max((U(index_U_max2,j,t) - MU(index_U_max2,j,t)),0) + K(index_U_max2,j,t);      %% Backlog Queue 2    (5) pag 481
                    else
                        U(i,j,t+1) = max(U(i,j,t),0) + K(i,j,t);           % Backlog Queue 2    (5) pag 481
                    end
                end
                
            end
        end                                                                % for each M-server
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%check points
    else
        bb=bb+1;
        UU1=zeros(N,M,TT);
        a_temp = zeros(1,M); 
        for jtest1=1:M
            if (a(jtest1)==1)
                UU1(:,jtest1,t)=U(:,jtest1,t);
            else
                UU1(:,jtest1,t)=nan;
            end
        end
        %if (t>1)
        %UU(~UU)=nan;
        %end
        for i= 1 : N
            
            if(W(i,t) > V.*alpha(i))   R(i,t) = 0;
                
            else
                R(i,t) = A(i,t);
            end
            
            %list U of active (ON) servers.
            
            [U2, indexu2]=min(UU1(i,:,t));% page 483 routing
            
            %[U2, indexu2]=min(U(i,:,t));                                   % page 483 routing
            if (W(i,t) > U2)
                K(i,indexu2,t) = a(indexu2)*W(i,t);
            else
                K(i,:,t)=0;
            end
            UU1(:,indexu2,t)=nan;
            %a(i,j) = (a(i,j) > 0.5);                                      % garantee to set more than half to be -on
            % fix number of active server, they are not M , they are
            % obviously less
            % if ((W(i,t)-sum(K(i,:,t)))>=0)
            
            % W(i,t+1) = W(i,t) - sum(a(i,j).*K(i,:,t)) + R(i,t);
            W(i,t+1) = W(i,t) - (a*K(i,:,t)') + R(i,t);                %% Backlog Queue 1    (2) pag 481 consider only active server
            % else
            %W(i,t+1)=R(i,t);
        end
        
        
        for j= 1 : M
            a_temp(j)=a(j);
            
             if (a1(j)==0)                                                  % for Dynamic server
                 
                 [U_max3, index_U_max3]=max(U(:,j,t));                           % 0=<index_U_max<=N
                 
                 for zzz=1:length(F3)
                     temp3(index_U_max3,j,zzz) = U_max3.*(min(mu3(index_U_max3,j,zzz),U_max3)) - V.*beta.*(P_cpu3(index_U_max3,j,zzz));                  % Resource allocation pag 483
                     
                 end
                 [value3, index3]=max(temp3(index_U_max3,j,:));                % 1=<index<=length(F)
                 
                 if (((min(mu3(index_U_max3,j,index3),U_max3))==0) && (P_cpu3(index_U_max3,j,index3)==0))
                     cc=cc+1;
                     a(j)=0;
                     MU(index_U_max3,j,t) = 0;
                     P(j,t)= 0;
                     if (t==1)
                         Switch(j,t)=k_e*(F3(index3))^2;
                         Freqtemp(j,2)=F3(index3);
                     else
                         Switch(j,t)=k_e*(F3(index3)-Freqtemp(j,t-1))^2;
                         Freqtemp(j,t)=F3(index3);
                     end
                 else
                     a(j)=1;
                     MU(index_U_max3,j,t) = min(mu3(index_U_max3,j,index3),U_max3);
                     P(j,t)= P_cpu3(index_U_max3,j,index3);
                     if (t==1)
                         Switch(j,t)=k_e*(F3(index3))^2;
                         Freqtemp(j,2)=F3(index3);
                     else
                         Switch(j,t)=k_e*(F3(index3)-Freqtemp(j,t-1))^2;
                         Freqtemp(j,t)=F3(index3);
                     end
                 end
            
                 for i= 1 : N
                     if (i==index_U_max3)
                         U(index_U_max3,j,t+1) = max(U(index_U_max3,j,t) - MU(index_U_max3,j,t),0) + K(index_U_max3,j,t);    % Backlog Queue 2    (5) pag 481
                     else
                         U(i,j,t+1) = max(U(i,j,t),0) + K(i,j,t);           % Backlog Queue 2    (5) pag 481
                     end
                 end
                 
            % for STATIC server % for STATIC server % for STATIC server % for STATIC server % for STATIC server % for STATIC server
            
             elseif (a1(j)==1)                                             % for STATIC server
                 
                [U_max4 , index_U_max4]=max(U(:,j,t));                     % 0=<index_U_max<=N
                MU(index_U_max4,j,t) = min(mu2(index_U_max4,j, length(F2)),U_max4);
                %P(j,t)= P_cpu2(index_U_max4,j,length(F2));
                P(j,t)= P_max2(index_U_max4,j,length(F2));
                a(j)=1;                                                     % STATIC SERVER IS ALWAYS ON
                if (t==1)
                    Switch(j,t)=k_e*(F2(length(F2)))^2;
                    Freqtemp(j,2)=F2(length(F2));
                else
                    Switch(j,t)=k_e*(F2(length(F2))-Freqtemp(j,t-1))^2;
                    Freqtemp(j,t)=F2(length(F2));
                end
                for i= 1 : N
                    if (i==index_U_max4)
                        U(index_U_max4,j,t+1) = max(U(index_U_max4,j,t) - MU(index_U_max4,j,t),0) + K(index_U_max4,j,t);    % Backlog Queue 2    (5) pag 481
                    else
                        U(i,j,t+1) = max(U(i,j,t),0) + K(i,j,t);           % Backlog Queue 2    (5) pag 481
                    end
                end
               
             end  % end of for Dynamic server
             
            % U(index_U_max,j,t+1) = max(U(index_U_max,j,t) - MU(index_U_max,j,t),0) + a(i,j).*K(index_U_max,j,t);      %% Backlog Queue 2    (5) pag 481
            % if
        end                                                                % for each M-server
        
        
        %*************************VM_MIGRATION******************
        on_servers = zeros(N+1,M);                                                 % on-server in each different frame T
        off_servers = zeros(N+1,M);                                                % off-server in each different frame T
        sleeping_servers= zeros(N+1,M);
        
        for jj = 1 : M
            for ii = 1 : N
                if ((a_temp(jj)==0) && (a(jj)==1)) || (((a_temp(jj)==1) && (a(jj)==1)))
                    on_servers(1,jj) = jj; %always ON servers + OFF--->> ON servers
                    on_servers(ii+1,jj) =  U(ii,jj,t);
                    %M1 = M1 + 1;
                elseif((a_temp(jj)==1) && (a(jj)==0))
                    off_servers(1,jj) = jj;
                    off_servers(ii+1,jj) = U(ii,jj,t);
                    %M2 = M2 + 1;
                    
                elseif ((a_temp(jj)==0) && (a(jj)==0))
                    sleeping_servers(1,jj)=jj;  %the servers which are never ON, ALWAYS OFF
                    sleeping_servers(ii+1,jj) =  U(ii,jj,t);
                    %M3=M3+1;
                end
                
            end
            
        end
                
        
        
%         on_servers = deletezero(on_servers);
%        [x1, x2]=size(on_servers);
%         i=0; j=0;
%         for j=1:x2
%             %for j=1:N
%             if (on_servers(2:1:(N+1),j)==0)
%               
%                 if ((a_temp(j)==1) && (a(j)==0) && (a1(j)==0))
%                     off_servers(1,j) = j;
%                     off_servers(2:1:(N+1),j) = 0;
%                 elseif ((a_temp(jj)==0) && (a(jj)==0) && (a1(j)==0))
%                     sleeping_servers(1,j)=j;  %the servers which are never ON, ALWAYS OFF
%                     sleeping_servers(2:1:(N+1),j) =  0;
%                 end
%             
%             end
%         end
        on_servers = deletezero(on_servers);
        on_servers = deletezero(on_servers);
        off_servers= deletezero(off_servers);
        sleeping_servers=deletezero(sleeping_servers);
        [x1, x2]=size(on_servers);
        M1=x2;
        [y1, y2]=size(off_servers);
        M2=y2;
        [z1, z2]=size(sleeping_servers);
        M3=z2;
        
        L_off(bb)=sum(sum(off_servers(2:N+1,:))) ;                                  %  the total migrated backlogs of all aplication in bits from the off servers
        for jk = 1 : M2
            for iii = 1 : N
                
                
                [a_tempON, indexON] = min(on_servers(iii+1,:));
                tempOnserver=off_servers(iii+1,jk);
                on_servers(iii+1,indexON) = on_servers(iii+1,indexON)+tempOnserver;
                off_servers(iii+1,jk)=0;
                
                
                
            end
        end
        
        
        a=zeros(1,M);
        for dom=1:M1
            if (on_servers(1,dom)>=1) && (a1(on_servers(1,dom))==0) % dynamic server
                a(on_servers(1,dom))=1;
                for i= 1 : N
                    U(i,on_servers(1,dom),t+1)=on_servers(i+1,dom);
                end
            elseif (on_servers(1,dom)>=1) && (a1(on_servers(1,dom))==1) % static server
                a(on_servers(1,dom))=1;
                for i= 1 : N
                    U(i,on_servers(1,dom),t+1)=on_servers(i+1,dom);
                end
            end
        end
        for dom=1:M2
            a(off_servers(1,dom))=0;
            for i= 1 : N
                U(i,off_servers(1,dom),t+1)=0;
            end
        end
                for dom=1:M3
                   a(sleeping_servers(1,dom))=0;
               end
        
        
        MM1(bb)=M1;
        MM2(bb)=M2;
        MM3(bb)=M3;
        M1=0;
        M2=0;
        M3=0;
    end % t=nT
    
    
end % biggggg for TT


for ii= 1:N
    r(ii) = mean(R(ii,:));
    ww(ii) = mean(W(ii,:));                                           % average occupation for each applications before router
    for jj = 1:M
    e(jj) =   mean(P(jj,:));                                      %pag482 (7)
    SW(jj)= mean (Switch(jj,:));
    AU(ii,jj) = mean(U(ii,jj,:));                                     % average occupation of service rate for application i in sererver j in block 2
    KU(ii,jj) = sum(K(ii,jj,:));                                     % average occupation of backlogs (in queue) of block 2 (queues in servers)
    end
    
end

KUU=0;

    for i= 1:N
        for j= 1:M
            KUU=KU(i,j)+KUU;
        end
    end

KUU=KUU/(M*N*TT);
for i= 1:N
Delay(i)=ww(i)/ r(i);  % delay takes for each application queue before routing

end


Avg_Delay1=mean(Delay(:)); % average delay of BLOCK 1 (before routing to the servers)



for i= 1:N
for j= 1:M
%Delay2(i,j)=AU(i,j)./ KU(i,j);  % delay takes for application i in server j in block 2
Delay2(i,j)=AU(i,j)./ KUU;  % delay takes for application i in server j in block 2
end
end



for j= 1:M
Avg_SDelay(j)=mean(Delay2(:,j)); % average delay takes for server j in block 2
end


Avg_Delay2=mean(Avg_SDelay(:));  % average dalay of BLOCK 2 (after routing to the servers and after servers )



TOT_DELAY=Avg_Delay1+Avg_Delay2;  %average total delay of WHOLE SYSTEM
UUTILITY = alpha(1).*sum(r(:)) - beta.*sum(e(:)); % UTILITY OF THE SYSTEM



%AVG_UTILITY = mean(UUTILITY(:));   % total utility of the SYSTEM





















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%for t= 1:TT-1
% 
% for ii= 1:N
%     r(ii) = mean(R(ii,:));                                             % throughput 4 each ii
%     ww(ii) = mean(W(ii,:));                                            % average occupation for each applications before router
%     
%     for jj = 1:M
%         e(jj) = mean(P(jj,:));                                      %pag482 (7)
%         AU(ii,jj) = mean(U(ii,jj,:));                                   % average occupation of backlogs applications of each server
%         avg_MU(ii,jj)= mean(MU(ii,jj,:));
%     end
%     
% end
% %end
% www=mean(ww(:));                                                           % average occpuation of all applications
% for i = 1 : N
%     AAU(i)=mean(AU(i,:));                                                       % average occupation of each backlogs application in the system
%     avg_AMU(i)=mean(avg_MU(i,:));
% end
% 
% for j = 1 : M
%     SAU(j)=mean(AU(:,j));                                                       % average occupation of total backlogs applications in each  server j
%     avg_SMU(jj)=mean(avg_MU(:,j));
% end
% 
% TAU=mean(AAU(:));                                                          % average occupation for : server for : applications
% Tavg_MU=mean(avg_AMU(:));
% 
% 
% AS = www + TAU ;                                                           % total average occupations of the system
% 
% Avg_DELAY =  AS/mean(r(:));                                                % average delay of the system
% 
% AAS(:) = ww(:) + AAU(:) ;                                                  % total average occupations of the system for each application
% 
% Avg_ADELAY(:) = AAS(:)./r(:);                                              % average delay of the system for each application
% 
% Avg_SDELAY(:) = SAU(:)./avg_SMU(:);
% 
% UTILITY = alpha(1)*sum(r(:)) - beta*sum(e(:));
% 
% 
































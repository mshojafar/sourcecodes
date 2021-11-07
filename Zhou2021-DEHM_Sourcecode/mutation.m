%为了保证多样性，在产生新的种群个体的过程中，产生的nrandI个互不相等的随机数，与i皆不相等；
%即：每产生的第 i 个新个体所用的随机选到的nrandI个旧个体不能是第 i 个旧个体。
function V=mutation(X,bestX,F,flag,finalx)
% function V=mutation(X,bestX,F,lamda,u,flag,finalx)
NP=size(X,1);
for i=1:NP
    %在[1 NP]中产生nrandI个互不相等的随机数，且与i皆不相等
    nrandI=4;
    r=randi([1,NP],1,nrandI);
    for j=1:nrandI
        equalr(j)=sum(r==r(j));
    end
    equali=sum(r==i);
    equalall=sum(equalr)+equali;
    while(equalall>nrandI) %若产生的随机数有相等的或与i相等的——需要重新生成随机数
        r=randi([1,NP],1,nrandI);
        for j=1:nrandI
            equalr(j)=sum(r==r(j));
        end
        equali=sum(r==i);
        equalall=sum(equalr)+equali;
    end
    if flag==0

         V(i,:)=X(i,:)+F.*(finalx-X(r(1),:))+F.*(bestX-X(r(2),:));   
    else
         temp=0.2*NP;
         k1=0;k2=0;
         while k1==k2
         k1=ceil(temp*rand);
         k2=ceil(temp*rand);
    end
      V(i,:)=X(i,:)+F.*(X(k1,:)-X(r(1),:))+F.*(X(k2,:)-X(r(2),:));   
    end
end

%为了保证多样性，在产生新的种群个体的过程中，产生的nrandI个互不相等的随机数，与i皆不相等；
%即：每产生的第 i 个新个体所用的随机选到的nrandI个旧个体不能是第 i 个旧个体。

function V=mutation(X,bestX,F,lamda,u)
NP=length(X);
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
    V(i,:)=u.*X(r(1),:)+(1-u).*bestX+F.*(lamda.*(X(r(2),:)-X(r(3),:))+(1-lamda).*(bestX-X(r(4),:)));   
end

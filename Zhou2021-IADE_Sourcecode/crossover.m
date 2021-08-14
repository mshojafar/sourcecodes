function [U,H]=crossover(X,V,CR)
[NP,Dim]=size(X);
for i=1:NP         %NP：种群个数；Dim:维度；V：变异向量；U：交叉向量；
    jRand=randi([1,Dim]);  %jRand∈[1,Dim]
    for j=1:Dim
        k=rand;
        if k<=CR||j==jRand  %j==jRand是为了确保至少有一个U(i,j)=V(i,j)
            U(i,j)=V(i,j);
        else
            U(i,j)=X(i,j);
        end
        if k>CR
            H(i,j)=X(i,j);
        else
            H(i,j)=V(i,j);
        end        
    end
end
end

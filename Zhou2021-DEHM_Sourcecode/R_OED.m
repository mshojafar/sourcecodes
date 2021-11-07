function [finalx]=R_OED(x,fitness,bestpop ,func_num,fhd)
[NP,d]=size(x);
level=2;
OO = OED(level,d)+1; %生成了xx个xx维的组合
[~,fx_sort] = sort(fitness); %升序
for i=1:level
    xx(i,:) = x(fx_sort(i),:);
end

x_test = zeros(size(OO,1),d);
% x_test：分类组合后的维度取出后存放
for i=1:size(OO,1)
    for j=1:d
        x_test(i,j) = xx(OO(i,j),j);  
    end
end

%xtemp：分类组合后的维度与背景矢量相结合形成新的位置
for i=1:size(x_test,1)
    xtemp(i,:)=x_test(i,:);
%     fx_test(i,:)=fobj(xtemp(i,:));
end
fx_test=(feval(fhd,xtemp',func_num))';
[~,ind]=min(fx_test);
xpop=xtemp(ind,:);
% 拼组合的过程
final_x = zeros(1,d);

for i = 1:d
    ssum = zeros(level,1);
    for jj = 1:level
        ssum(jj,1) = sum(fx_test(find(OO(:,i)==jj)));
    end
    [a,b]=min(ssum);
    final_x(1,i)=xx(b(1),i);
end

for k=1:d
    finalx(:,k)=final_x(:,k);
end
% ffinal_x = fobj(finalx);
ffinal_x=feval(fhd,finalx',func_num);

%ffinal_x < Fbest
%fes=fes+1;
if ffinal_x>min(fx_test)
    finalx=xpop;
    %   min(fx_test)<Fbest
end
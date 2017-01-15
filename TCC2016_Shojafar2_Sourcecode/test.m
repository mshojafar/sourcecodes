M=46;
Q=4;
TT1=zeros(M,Q+1);
TT_Q=zeros(1,Q+1);
t_opt_wokloads1=zeros(1,Q+1);

for i=1:M
    for j=1:Q
        TT1(i,j)=mean(Total_Time(:,i,j));
    end
end
for j=1:Q+1
    %TT_Q(i)=mean(TT1(:,j));
    t_opt_wokloads1(j)=mean(t_opt_wokloads(:,j));
end
t_opt_wokloads(~any(isnan(t_opt_wokloads)),:)
t_opt_wokloads(isnan(t_opt_wokloads)) = [];
bar(t_opt_wokloads1(:)/sum(t_opt_wokloads1(:)));
$\overline{L_{tot}}$=593.3, $max(L_{tot})$=1667, M=50, T=3,T_t=5

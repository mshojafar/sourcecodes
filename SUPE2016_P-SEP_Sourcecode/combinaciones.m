clc;
n=3;
m=2;
j=1;
h=0;
for i=1:1:(n^m)
    if (i>0)&&(i<7)
        if (h>1)&&(h<4)
            j=j+1;
            l=[h,j];
            display(l);
        end
        l=[i,j];
        h=h+1;
        display(l);
    end
end

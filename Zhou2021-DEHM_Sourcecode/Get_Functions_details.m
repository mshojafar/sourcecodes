function [lb,ub,fobj] = Get_Functions_details(F)



switch F
    case 1
        fobj = @F1;
        lb=-100;
        ub=100;
        
        
    case 2
        fobj = @F2;
        lb=-100;
        ub=100;
        
        
    case 3
        fobj = @F3;
        lb=-100;
        ub=100;
        
        
    case 4
        fobj = @F4;
        lb=-100;
        ub=100;
        
        
    case 5
        fobj = @F5;
        lb=-100;
        ub=100;
        
        
    case 6
        fobj = @F6;
        lb=-100;
        ub=100;
        
        
    case 7
        fobj = @F7;
        lb=-100;
        ub=100;
        
        
    case 8
        fobj = @F8;
        lb=-100;
        ub=100;
        
        
    case 9
        fobj = @F9;
        lb=-100;
        ub=100;
        
        
    case 10
        fobj = @F10;
        lb=-32;
        ub=32;
        
        
    case 11
        fobj = @F11;
        lb=-100;
        ub=100;
        
        
    case 12
        fobj = @F12;
        lb=-50;
        ub=50;
        
        
    case 13
        fobj = @F13;
        lb=-100;
        ub=100;
        
        
    case 14
        fobj = @F14;
        lb=-100;
        ub=100;
        
        
    case 15
        fobj = @F15;
        lb=-100;
        ub=100;
        
        
    case 16
        fobj = @F16;
        lb=-100;
        ub=100;
        
        
    case 17
        fobj = @F17;
        lb=-100;
        ub=100;
        
        
    case 18
        fobj = @F18;
        lb=-100;
        ub=100;
        
        
    case 19
        fobj = @F19;
        lb=-100;
        ub=100;
        
    case 20
        fobj = @F20;
        lb=-100;
        ub=100;
        
    case 21
        fobj = @F21;
        lb=-100;
        ub=100;
        
        
    case 22
        fobj = @F22;
        lb=-100;
        ub=100;
        
        
    case 23
        fobj = @F23;
        lb=-100;
        ub=100;
        
    case 24
        fobj = @F24;
        lb=-100;
        ub=100;
        
    case 25
        fobj = @F25;
        lb=-100;
        ub=100;
        
    case 26
        fobj = @F26;
        lb=-100;
        ub=100;
        
    case 27
        fobj = @F27;
        lb=-100;
        ub=100;
        
    case 28
        fobj = @F28;
        lb=-100;
        ub=100;
        
    case 29
        fobj = @F29;
        lb=-100;
        ub=100;
        
    case 30
        fobj = @F30;
        lb=-100;
        ub=100;
        
end
end

% F1

function oo = F1(x)
oo=sum(x.^2);
end

% F2

function oo = F2(x)
oo=sum(abs(x))+prod(abs(x));
end

% F3

function oo = F3(x)
dim=size(x,2);
oo=0;
for i=1:dim
    oo=oo+sum(x(1:i))^2;
end
end

% F4

function oo = F4(x)
oo=max(abs(x));
end

% F5

function oo = F5(x)
dim=size(x,2);
oo=sum(100*(x(2:dim)-(x(1:dim-1).^2)).^2+(x(1:dim-1)-1).^2);
end

% F6

function oo = F6(x)
oo=sum(abs((x+.5)).^2);
end

% F7

function oo = F7(x)
dim=size(x,2);
oo=sum([1:dim].*(x.^4))+rand;
end

% F8

function oo = F8(x)
oo=abs(sum(-x.*sin(sqrt(abs(x)))));
end

% F9

function oo = F9(x)
dim=size(x,2);
oo=sum(x.^2-10*cos(2*pi.*x))+10*dim;
end

% F10

function oo = F10(x)
dim=size(x,2);
oo=-20*exp(-.2*sqrt(sum(x.^2)/dim))-exp(sum(cos(2*pi.*x))/dim)+20+exp(1);
end

% F11

function oo = F11(x)
dim=size(x,2);
oo=sum(x.^2)/4000-prod(cos(x./sqrt([1:dim])))+1;
end

% F12

function oo = F12(x)
dim=size(x,2);
oo=(pi/dim)*(10*((sin(pi*(1+(x(1)+1)/4)))^2)+sum((((x(1:dim-1)+1)./4).^2).*...
    (1+10.*((sin(pi.*(1+(x(2:dim)+1)./4)))).^2))+((x(dim)+1)/4)^2)+sum(Ufun(x,10,100,4));
end

% F13

function oo = F13(x)
dim=size(x,2);
oo=.1*((sin(3*pi*x(1)))^2+sum((x(1:dim-1)-1).^2.*(1+(sin(3.*pi.*x(2:dim))).^2))+...
    ((x(dim)-1)^2)*(1+(sin(2*pi*x(dim)))^2))+sum(Ufun(x,5,100,4));
end

% F14

function oo = F14(x)
global o M_10 M_30 M_50 M_100
D=size(x,2);
ok = o(1:D);
x = x - repmat(ok,size(x,1),1);
f=sum(x.^2,2);
oo=20-20.*exp(-0.2.*sqrt(f./D))-exp(sum(cos(2.*pi.*x),2)./D)+exp(1);
end

% F15

function oo = F15(x)
global o M_10 M_30 M_50 M_100
D=size(x,2);
ok = o(1:D);
x = x - repmat(ok,size(x,1),1);
oo=max(abs(x),[],2);
end

% F16

function oo = F16(x)
global o M_10 M_30 M_50 M_100
D=size(x,2);
ok = o(1:D);
x = x - repmat(ok,size(x,1),1);
oo=sum(abs(x),2);
end

% F17

function oo = F17(x)
global o M_10 M_30 M_50 M_100
D=size(x,2);
ok = o(1:D);
x = x - repmat(ok,size(x,1),1);
f=1;
for i=1:D
    f=f.*(cos(x(:,i)/sqrt(i)));
end
oo=sum(x.^2,2)/4000-f+1;
end

% F18

function oo = F18(x)
global o M_10 M_30 M_50 M_100
D=size(x,2);
ok = o(1:D);
x = x - repmat(ok,size(x,1),1);
f=1;
for i=1:D
    f=f.*(cos(x(:,i)/sqrt(i)));
end
oo=sum(x.^2,2)/4000-f+1;
end

% F19

function oo = F19(x)
global o M_10 M_30 M_50 M_100
D=size(x,2);
ok = o(1:D);
x = x - repmat(ok,size(x,1),1);
oo = sum(abs(x).^0.5 + 2*sin(x.^3), 2);
end

% F20

function oo = F20(x)
global o M_10 M_30 M_50 M_100
D=size(x,2);
ok = o(1:D);
x = x - repmat(ok,size(x,1),1);
sum1 =  sum(0.5 + (sin((x(:,1:D-1).^2 + x(:,2:D).^2).^0.5).^2 - 0.5)./(1+0.001*(x(:,1:D-1).^2 + x(:,2:D).^2).^0.5).^2,2);
sum2 =  0.5 + (sin((x(:,D).^2 + x(:,1).^2).^0.5).^2 - 0.5)./(1+0.001*(x(:,D).^2 + x(:,1).^2).^0.5).^2;
oo = sum1 + sum2;
end

% F21

function oo = F21(x)
global o M_10 M_30 M_50 M_100
D=size(x,2);
ok = o(1:D);
if(D==10)
    M = M_10;
end
if(D==30)
    M = M_30;
end
if(D==50)
    M = M_50;
end
if(D==100)
    M = M_100;
end

x = x - repmat(ok,size(x,1),1);
x = (M*x')';
oo = sum(x.^2-10.*cos(2.*pi.*x)+10,2);
end

% F22

function oo = F22(x)
global o M_10 M_30 M_50 M_100
D=size(x,2);
ok = o(1:D);
if(D==10)
    M = M_10;
end
if(D==30)
    M = M_30;
end
if(D==50)
    M = M_50;
end
if(D==100)
    M = M_100;
end

x = x - repmat(ok,size(x,1),1);
x = (M*x')';
oo = sum(100.*(x(:,1:D-1).^2-x(:,2:D)).^2+(x(:,1:D-1)-1).^2,2);
end

% F23

function oo = F23(x)
global o M_10 M_30 M_50 M_100
D=size(x,2);
ok = o(1:D);
if(D==10)
    M = M_10;
end
if(D==30)
    M = M_30;
end
if(D==50)
    M = M_50;
end
if(D==100)
    M = M_100;
end
x = x - repmat(ok,size(x,1),1);
x = (M*x')';
f=sum(x.^2,2);
oo=20-20.*exp(-0.2.*sqrt(f./D))-exp(sum(cos(2.*pi.*x),2)./D)+exp(1);
end

function oo=Ufun(x,a,k,m)
oo=k.*((x-a).^m).*(x>a)+k.*((-x-a).^m).*(x<(-a));
end

%F24
function oo = F24(x)
global o M_10 M_30 M_50 M_100
D=size(x,2);
ok = o(1:D);
if(D==10)
    M = M_10;
end
if(D==30)
    M = M_30;
end
if(D==50)
    M = M_50;
end
if(D==100)
    M = M_100;
end
x = x - repmat(ok,size(x,1),1);
x = (M*x')';
oo=max(abs(x),[],2);
end

%F25
function oo = F25(x)
global o M_10 M_30 M_50 M_100
D=size(x,2);
ok = o(1:D);
if(D==10)
    M = M_10;
end
if(D==30)
    M = M_30;
end
if(D==50)
    M = M_50;
end
if(D==100)
    M = M_100;
end
x = x - repmat(ok,size(x,1),1);
x = (M*x')';
sum1 = sum(abs(x),2);
oo= sum1;
end

%F26
function oo = F26(x)
global o M_10 M_30 M_50 M_100
D=size(x,2);
ok = o(1:D);
if(D==10)
    M = M_10;
end
if(D==30)
    M = M_30;
end
if(D==50)
    M = M_50;
end
if(D==100)
    M = M_100;
end

x = x - repmat(ok,size(x,1),1);
x = (M*x')';
f=1;
for i=1:D
    f=f.*(cos(x(:,i)/sqrt(i)));
end
oo=sum(x.^2,2)/4000-f+1;
end

%F27
function oo = F27(x)
global o M_10 M_30 M_50 M_100
D=size(x,2);
ok = o(1:D);
if(D==10)
    M = M_10;
end
if(D==30)
    M = M_30;
end
if(D==50)
    M = M_50;
end
if(D==100)
    M = M_100;
end
x = x - repmat(ok,size(x,1),1);
x = (M*x')';
g1 = -sum(abs(x),2) + 1;
g2 = sum(x.^2, 2) - 100*D;
sum1 = sum(100.*(x(1,1:D-1).^2-x(1,2:D)).^2);
multi = 1;
for i=1:D
    multi = multi.* sin((x(:,i)-1)*pi).^2;
end
h = sum1 + multi;
x=(abs(x)<0.5).*x+(abs(x)>=0.5).*(round(x.*2)./2);
oo=sum(x.^2-10.*cos(2.*pi.*x)+10,2);
end

%F28
function oo = F28(x)
global o M_10 M_30 M_50 M_100
D=size(x,2);
ok = o(1:D);
if(D==10)
    M = M_10;
end
if(D==30)
    M = M_30;
end
if(D==50)
    M = M_50;
end
if(D==100)
    M = M_100;
end
x = x - repmat(ok,size(x,1),1);
x = (M*x')';
oo = sum(abs(x).^0.5 + 2*sin(x.^3), 2);
end


%F29
function oo = F29(x)
global o M_10 M_30 M_50 M_100
D=size(x,2);
ok = o(1:D);
if(D==10)
    M = M_10;
end
if(D==30)
    M = M_30;
end
if(D==50)
    M = M_50;
end
if(D==100)
    M = M_100;
end
x = x - repmat(ok,size(x,1),1);
x = (M*x')';
sum1 =  sum(0.5 + (sin((x(:,1:D-1).^2 + x(:,2:D).^2).^0.5).^2 - 0.5)./(1+0.001*(x(:,1:D-1).^2 + x(:,2:D).^2).^0.5).^2,2);
sum2 =  0.5 + (sin((x(:,D).^2 + x(:,1).^2).^0.5).^2 - 0.5)./(1+0.001*(x(:,D).^2 + x(:,1).^2).^0.5).^2;
oo = sum1 + sum2;
end


%F30
function oo = F30(x)
global o M_10 M_30 M_50 M_100
D=size(x,2);
ok = o(1:D);
if(D==10)
    M = M_10;
end
if(D==30)
    M = M_30;
end
if(D==50)
    M = M_50;
end
if(D==100)
    M = M_100;
end
x = x - repmat(ok,size(x,1),1);
x = (M*x')';
sum1 =  sum(0.5 + (sin((x(:,1:D-1).^2 + x(:,2:D).^2).^0.5).^2 - 0.5)./(1+0.001*(x(:,1:D-1).^2 + x(:,2:D).^2).^0.5).^2,2);
sum2 =  0.5 + (sin((x(:,D).^2 + x(:,1).^2).^0.5).^2 - 0.5)./(1+0.001*(x(:,D).^2 + x(:,1).^2).^0.5).^2;
oo = sum1 + sum2.^2;
end
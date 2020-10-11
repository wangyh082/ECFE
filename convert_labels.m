function labels=convert_labels(g)
n=length(g);
x=[1:n];
W=sparse([x,g],[g,x],ones(1,2*n));
[r,c,v]=find(W);
r=r';
c=c';
temp=diff(c);
jin=find([1,temp,1]);

count=1;
labels=zeros(1,n);
v=zeros(1,n);

s=zeros(1,n);
sp=1;

for i=1:n
    if v(i)==1
        continue;
    end
    sp=1;
    s(sp)=i;
    sp=sp+1;
    while sp~=1
        sp=sp-1;
        t=s(sp);
        labels(t)=count;
        v(t)=1;
        for j=jin(t):jin(t+1)-1
            k=r(j);
            if v(k)~=1
                s(sp)=k;
                sp=sp+1;
            end
        end
    end
    count=count+1;
end


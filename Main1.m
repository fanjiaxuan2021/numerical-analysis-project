%给定p,q,f,x00,x01,x10,x11,a,b,ui,ui0,tol
k=16;
C=4;
test=1;
f=@(x) 1000;
[gl,gr,g0l,g0r] = green(a,b,x00,x01,x10,x11);
s=gl(a)*g0r(a)-gr(a)*g0l(a);
phil=@(x) (p(x)*g0r(x)+q(x)*gr(x))/s;
phir=@(x) (p(x)*g0l(x)+q(x)*gl(x))/s;
new_node=Node.empty(1000,0); %1000待定
new_node(1)=Node(a,b);
new_node(1).Value=localinnerproduct(a,b,gl,gr,phil,phir,f,k);
new_node(1).Lambda=[0,0];
n=1; %所有node个数
r=1;
new_node(1).Monitor=monitor(a,b,phil,phir,gl,gr,f,k); %monitor的定义还需进一步完善以及'f'是什么
Sdiv=new_node(1).Monitor;
flag0=1;
flag1=1;

while flag0==1 || flag1==1
    %以下是细分算法未考虑list&重复环节
    flag0=0;
    flag1=0;
    list=[];
    for s=1:n
        if isempty(new_node(s).Location)
                list=[list,s];
        end
    end
    m=length(list);
    for i=1:m
        new_node(list(i)).Value=localinnerproduct(new_node(list(i)).Name(1),new_node(list(i)).Name(2),gl,gr,phil,phir,f,k);
    end
    %upward sweep
    for i=n:1
        if ~isempty(new_node(i).Location)
            X=new_node(new_node(i).Location(1)).Value;
            Y=new_node(new_node(i).Location(1)+1).Value;
            new_node(i).Value=upward(X,Y);
        end
    end
    %downward sweep
    for i=1:n
        if ~isempty(new_node(i).Location)
            X=new_node(new_node(i).Location(1)).Value;
            Y=new_node(new_node(i).Location(1)+1).Value;
            T=new_node(i).Lambda;
            [new_node(new_node(i).Location(1)).Lambda,new_node(new_node(i).Location(1)+1).Lambda]=downward(X,Y,T);
        end
    end
    %给leafnode排序
    Q=zeros(1,m);
    for i=1:m
        t=1;
        for j=1:m
            if new_node(list(i)).Name(1)>new_node(list(j)).Name(1)
                t=t+1;
            end
        end
        Q(t)=list(i);
    end
    %sigma的近似在节点的值
    point=zeros(1,m*k);
    sigma=zeros(1,m*k);
    for i=1:m
        a=new_node(list(i)).Name(1);b=new_node(list(i)).Name(2);
        tau=cheb(a,b,k);
        P=local(a,b,phil,phir,gl,gr,k);
        pl=eva(phil,a,b,k);
        pr=eva(phir,a,b,k);
        f0=eva(f,a,b,k);
        v1=P\f0;
        v2=P\pl;
        v3=P\pr;
        for j=1:k
            point(i*k-k+j)=tau(j);
            sigma(i*k-k+j)=v1(j)+new_node(list(i)).Lambda(1)*v2(j)+new_node(list(i)).Lambda(2)*v3(j); 
        end
    end
    %积分近似
    jl=zeros(1,m);jr=zeros(1,m);
    jl(1)=0;jr(m)=0;
    if m>=2
        for i=1:m-1
            jl(i+1)=jl(i)+new_node(Q(i)).Value(5)+new_node(Q(i)).Lambda(1)*new_node(Q(i)).Value(1)+new_node(Q(i)).Lambda(2)*new_node(Q(i)).Value(3);
        end
        for i=m:-1:2
            jr(i-1)=jr(i)+new_node(Q(i)).Value(6)+new_node(Q(i)).Lambda(1)*new_node(Q(i)).Value(2)+new_node(Q(i)).Lambda(2)*new_node(Q(i)).Value(4);
        end
    end
    %解在节点的近似
    u=zeros(1,m*k);
    u0=zeros(1,m*k);
    for i=1:m
        a=new_node(list(i)).Name(1);b=new_node(list(i)).Name(2);
        z=zeros(1,k);
        hl=eva(gl,a,b,k);
        hr=eva(gr,a,b,k);
        for j=1:k
            z(j)=sigma(i*k-k+j);
        end
        y1=il(k)*(hl.*z);y2=ir(k)*(hr.*z);
        for j=1:k
            u(i*k-k+j)=ui(p(i*k-k+j))+gr(p(i*k-k+j))*(jl(i)+y1(j))/s+gl(p(i*k-k+j))*(jr(i)+y2(j))/s;
            u0(i*k-k+j)=ui0(p(i*k-k+j))+g0r(p(i*k-k+j))*(jl(i)+y1(j))/s+g0l(p(i*k-k+j))*(jr(i)+y2(j))/s;
        end
    end

    S=zeros(1,m);
    for i=1:m
        S(i)=monitor(new_node(list(i)).Name(1),new_node(list(i)).Name(2),phil,phir,gl,gr,f,k);
    end
    Sdiv=max(max(S));
    %算test
    if r>1
        L1=zeros(1,m*k);L2=zeros(1,m*k);
        for i=1:m*k
            L1=abs(ff(point(i))-u(i));
            L2=abs(ff(point(i))+u(i));
        end
        test=max(max(L1))/max(max(L2));
    end
    if test>tol
        for i=1:m
            if new_node(list(i)).Monitor>=Sdiv
                l=new_node(list(i)).Name(1);
                r=new_node(list(i)).Name(1);
                new_node(n+1)=Node(l,(l+r)/2);
                new_node(n+2)=Node((l+r)/2,r);
                new_node(n+1).Parent=list(i);
                new_node(n+2).Parent=list(i);
                new_node(list(i)).Location=n+1;
                n=n+2;
                flag0=1;
            end
        end
    
        %以下是合并算法
        dp=[]; %被删父节点
        dc=[]; %被删子节点
        for j=1:m/2
            if new_node(list(2*j-1)).Monitor+new_node(list(2*j)).Monitor<(Sdiv/(2^16)) %16次方待定
                flag1=1;
                dp=[dp,2*j-1,2*j];
                if isempty(new_node(2*j-1).Location)
                    dc=[dc,new_node(2*j-1).Location,new_node(2*j-1).Location+1];
                end
                if isempty(new_node(2*j).Location)
                    dc=[dc,new_node(2*j).Location,new_node(2*j).Location+1];
                end
            end
            dq=[dp,dc];

            for u=1:(length(dq)/2-1) %计算所有变化子节点
                new_node(new_node(dq(2*u-1)).Parent).Location=[];
                for p=new_node(dq(2*u-1)).Parent+1:new_node(dq(2*u+1)).Parent-1
                    if isempty(new_node(p).Location)==0
                        new_node(p).Location=new_node(p).Location-2*u;
                    end
                end
            end
            for p=new_node(dq(end)).Parent+1:new_node(list(end)).Parent
                new_node(p).Location=new_node(p).Location-length(dq);
            end

            for v=1:(length(dp)-1) %计算所有变化父节点
                for p=new_node(dp(v)).Location+1:new_node(dp(v+1)).Location-1
                    new_node(p).Parent=new_node(p).Location-v;
                end
            end
            for p=new_node(dq(end)).Location+1:new_node(list(end)).Location+1
                new_node(p).Parent=new_node(p).Parent-length(dq);
            end

            for w=1:length(dq)
                new_node(dq(w))=[];
            end

        end
    end
    %定义用于算范数的函数ff
    c=zeros(1,m+1);
    for i=1:m
        c(i)=new_node(Q(i)).Name(1);
    end
    c(m+1)=b;
    ff=@(x) stepcheb(c,u,k,x);
    r=r+1;
end
for s=1:n
    list=[];
    if isempty(new_node(s).Location)
                list=[list,s];
    end
end
m=length(list);
for i=1:m
    new_node.Value=localinnerproduct(new_node(i).Name(1),new_node(i).Name(2),gl,gr,phil,phir,f,k);
end
for i=n:1
    if ~isempty(new_node(i).Location)
        X=new_node(new_node(i).Location(1)).Value;
        Y=new_node(new_node(i).Location(2)).Value;
        new_node(i).value=upward(X,Y);
    end
end
for i=1:n
    if ~isempty(new_node(i).Location)
        X=new_node(new_node(i).Location(1)).Value;
        Y=new_node(new_node(i).Location(2)).Value;
        T=new_node(i).Lamda;
        [new_node(new_node(i).Location(1)).Lambda,new_node(new_node(i).Location(2)).Lambda]=downward(X,Y,T);
    end
end
Q=zeros(1,m);
for i=1:m
    t=1;
    for j=1:m
        if new_node(i).Name(1)>new_node(j).Name(1)
            t=t+1;
        end
    end
    Q(t)=list(i);
end
point=zeros(1,m*k);
sigma=zeros(1,m*k);
for i=1:m
    a=new_node(list(i)).Name(1);b=new_node(list(i)).Name(2);
    tau=cheb(a,b,k);
    P=local(a,b,phil,phir,gl,gr,k);
    pl=eva(phil,a,b,k);
    pr=eva(phir,a,b,k);
    f0=eva(f,a,b,k);
    v1=P\f0;
    v2=P\pl;
    v3=P\pr;
    for j=1:k
        point(i*k-k+j)=tau(j);
        sigma(i*k-k+j)=v1(j)+new_node(list(i)).Lambda(1)*v2(j)+new_node(list(i)).Lambda(2)*v3(j);
    end
end
jl=zeros(1,m);jr=zeros(1,m);
jl(1)=0;jr(m)=0;
for i=1:m-1
    jl(i+1)=jl(i)+new_node(Q(i)).Value(5)+new_node(Q(i)).Lambda(1)*new_node(Q(i)).Value(1)+new_node(Q(i)).Lambda(2)*new_node(Q(i)).Value(3);
    jr(m-i)=jr(m-i+1)+new_node(Q(m-i+1)).Value(6)+new_node(Q(m-i+1)).Lambda(1)*new_node(Q(m-i+1)).Value(2)+new_node(Q(m-i+1)).Lambda(2)*new_node(Q(m-i+1)).Value(4);
end
u=zeros(1,m*k);
u0=zeros(1,m*k);
for i=1:m
    a=new_node(list(i)).Name(1);b=new_node(list(i)).Name(2);
    z=zeros(1,k);
    hl=eva(gl,a,b,k);
    hr=eva(gr,a,b,k);
    for j=1:k
        z(j)=sigma(i*k-k+j);
    end
    y1=il(k)*(hl.*z);y2=ir(k)*(hr.*z);
    for j=1:k
        u(i*k-k+j)=ui(p(i*k-k+j))+gr(p(i*k-k+j))*(jl(i)+y1(j))/s+gl(p(i*k-k+j))*(jr(i)+y2(j))/s;
        u0(i*k-k+j)=ui0(p(i*k-k+j))+g0r(p(i*k-k+j))*(jl(i)+y1(j))/s+g0l(p(i*k-k+j))*(jr(i)+y2(j))/s;
    end
end
plot(point,u)
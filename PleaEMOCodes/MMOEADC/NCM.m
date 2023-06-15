function C=NCM(X,r)
    [n,dim]=size(X);
    S=(1:n);
    D=zeros(n,n);
    for i=1:dim
        D=D+(abs(X(:,i)-repmat(X(:,i)',n,1))<=r(i));
    end

    K=0;
    while size(S)>0
        K=K+1;
        Q=[];
        C(K).p=[];
        s=randperm(length(S));
        x=S(s(1));
        Q=[Q x];
        C(K).p=[C(K).p x];
        while size(Q)>0
            ss=randperm(length(Q));
            y=Q(ss(1));
            B=find(D(y,:)==dim);
            T=setdiff(B,C(K).p);
            Q=[Q T];
            C(K).p=[C(K).p T];
            Q(ss(1))=[];
        end
        S=setdiff(S,C(K).p);
    end
end
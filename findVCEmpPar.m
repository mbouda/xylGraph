function  VC=findVCEmpPar(G,sp,p0,p1,nT,res,stPt,nEmbInit)

    dP=p0:0.5:p1;
    nP=size(dP,2);
    pe=zeros(nT,nP);
    parfor i=1:nP
        pel=zeros(nT,1);
        dpl=dP(i);
        for k=1:nT
            pel(k)=emboliseGraphEmp(G,sp,dpl,nEmbInit);
        end
        pe(:,i)=pel;
    end

    xp1=reshape(repmat(dP,[nT 1]),[nP*nT 1]);
    yp1=reshape(pe,[nP*nT 1]);

    [f,~]=fit(xp1,yp1,'1./(1 + exp(a*(x-b)))','StartPoint',stPt); %start point for sigma ASP units: [3 0]; otherwise [1 3]
    p05=fzero(@(x)f(x)-0.05,f.b);
    p95=fzero(@(x)f(x)-0.95,f.b);

    dP2=p05:(p95-p05)/(res-1):p95;
    pe2=zeros(nT,res);
    parfor i=1:res
        pel=zeros(nT,1);
        dpl=dP2(i);
        for k=1:nT
            pel(k)=emboliseGraphEmp(G,sp,dpl,nEmbInit);
        end
        pe2(:,i)=pel;
    end

    xp2=reshape(repmat(dP2,[nT 1]),[res*nT 1]);
    yp2=reshape(pe2,[res*nT 1]);

    %[f2,g2]=fit(xp2,yp2,'1./(1 + exp(a*(x-b)))','StartPoint',[1 3]);

    [f3,g3]=fit(cat(1,xp1,xp2),cat(1,yp1,yp2),'1./(1 + exp(a*(x-b)))','StartPoint',stPt); %start point for sigma ASP units: [3 0]; otherwise [1 3]


    VC.dP=cat(1,xp1,xp2);
    VC.pE=cat(1,yp1,yp2);
    VC.f=f3;
    VC.g=g3;
    VC.p12=fzero(@(x)f3(x)-0.12,f3.b);
    VC.p50=fzero(@(x)f3(x)-0.50,f3.b);
    VC.p88=fzero(@(x)f3(x)-0.88,f3.b);
end
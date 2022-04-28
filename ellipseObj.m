function J=ellipseObj(rb,targP88,ra,crds,th,r)
    
    incl=r<=(ra*rb)./sqrt((ra*sin(th)).^2+(rb*cos(th)).^2);
    
    crdG=crds(incl,:);

    D=distFun(crdG(:,1),crdG(:,2),size(crdG,1));
    conn=abs(D-1)<1e-6;
    [ei,ej]=find(triu(conn));
    G=graph(ei,ej);

    X=(-4:0.01:4)';
    F=normcdf(X,0,1);
    sp=spaps(F,X,0);

    p0=-3.5;
    p1= 3.5;
    nT=5000;
    res=64;
    nEmbInit=1;  %add initial embolism here if needed
    stPt=[0 3];
    
    VC=findVCEmpPar(G,sp,p0,p1,nT,res,stPt,nEmbInit);
    
    if VC.p88>5
        J=NaN;
    else
        J=VC.p88-targP88;
    end
        
end
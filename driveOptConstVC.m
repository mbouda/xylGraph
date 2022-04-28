function el=driveOptConstVC(nCores,ra,lb,ub,targP88,nV,nL)

    myPool=parpool(nCores);
    
    crds=honeycomb(nV,nL);
    [th,r]=cart2pol(crds(:,1),crds(:,2));
    
    rb=fzero(@(rb)ellipseObj(rb,targP88,ra,crds,th,r),[lb ub]);
    incl=r<=(ra*rb)./sqrt((ra*sin(th)).^2+(rb*cos(th) ).^2);

    crdG=crds(incl,:);
    D=distFun(crdG(:,1),crdG(:,2),size(crdG,1));
    conn=abs(D-1)<1e-6;
    [ei,ej]=find(triu(conn));
    G=graph(ei,ej);
    
    el.ra=ra;
    el.rb=rb;
    el.crds=crdsG;
    el.G=G;

    delete(myPool);
end
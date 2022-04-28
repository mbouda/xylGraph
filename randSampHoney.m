function randHG=randSampHoney(nVTarg,nG,nCores,nPts,nLay)

    rngs=rng('shuffle');
        %should save rngs for replicability

    crds=honeycomb(nPts,nLay);
    %have circle of radius=2140 filled with honeycomb points

    randHG=struct('NV',nVTarg,'nG',nG,'rngs',rngs,...
        'GS',struct('crds',cell(nG,1),'E',cell(nG,1),'G',cell(nG,1)));
    
    myPool=parpool(nCores);
    parfor i=1:nG
        inds=(1:nPts)';
        inShape=[];
        sel=1;
        inds=setdiff(inds,sel);
        selCrd=crds(sel,:);
        remCrd=crds(inds,:);

        d=sqrt((remCrd(:,1)-selCrd(1)).^2+(remCrd(:,2)-selCrd(2)).^2);
        el=inds(abs(d-1)<1e-6);

        for j=1:nVTarg-1
            inShape=cat(1,inShape,sel);
            sel=randsample(el,1);
            inds=setdiff(inds,sel);
            el=setdiff(el,sel);

            selCrd=crds(sel,:);
            remCrd=crds(inds,:);
            d=sqrt((remCrd(:,1)-selCrd(1)).^2+(remCrd(:,2)-selCrd(2)).^2);
            el=union(el,inds(abs(d-1)<1e-6));
        end
        inShape=cat(1,inShape,sel);
        %it's about 3 seconds per 100 nodes
        %now it's about 75 seconds per 100 nodes
        %with bigger crds array

        shapeCrds=crds(inShape,:);
        D=distFun(shapeCrds(:,1),shapeCrds(:,2),nVTarg); %can't do distFun like this with so many points
        conn=abs(D-1)<1e-6;
        [ei,ej]=find(triu(conn));
        E=[ei,ej];
        G=graph(ei,ej);
        
        GS(i).crds=shapeCrds;
        GS(i).G=G;
        GS(i).E=E;
    end
    randHG.GS=GS;
    delete(myPool);
end
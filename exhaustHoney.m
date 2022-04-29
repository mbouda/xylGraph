function exhGraphs=exhaustHoney(nPts,nLay,targSize)


    crds=honeycomb(nPts,nLay);
    D=distFun(crds(:,1),crds(:,2),nPts);

    conn=abs(D-1)<1e-6;
    [ei,ej]=find(triu(conn));
    G=graph(ei,ej);

    doneSet=cell(0,1);
    doneWgt=zeros(0,1);
    nodeSet={[1 2]};
    wghtSet=1;
    nVG=2;
    while nVG<targSize
        szNS=size(nodeSet);
        subSet=cell(szNS);
        dupSet=cell(szNS);
        for i=1:szNS(1)
            nodes=nodeSet{i};
            adj=setdiff(union(ej(ismember(ei,nodes)),ei(ismember(ej,nodes))),nodes);
            nAdj=size(adj,1);
            propSets=cat(2,repmat(nodes,[nAdj 1]),adj);

            match=idSimilar(propSets,crds,nAdj);
            dups=sum(match,2);
            subSet{i}=sort(propSets(~any(match,1)',:),2);
            dupSet{i}=dups(~any(match,1)');
        end

        for i=1:szNS(1)-1
            for j=i+1:szNS
                %pairwise comparisons between subSets
                %eliminate duplictaes
                [dup,dupii]=ismember(subSet{j},subSet{i},'rows');
                dupSet{i}(dupii(dup))=dupSet{i}(dupii(dup))+dupSet{j}(dup)+1;
                dupSet{j}=dupSet{j}(~dup);
                subSet{j}=subSet{j}(~dup,:);
                match=idSimSubs(subSet{i},subSet{j},crds);
                [k,l]=find(match);
                dupSet{i}(k)=dupSet{i}(k)+dupSet{j}(l)+1;
                dupSet{j}(l)=[];
                subSet{j}(any(match)',:)=[];
            end
        end

        for i=1:szNS(1)
            dupSet{i}=dupSet{i}+1;
            if ~any(dupSet{i})
                dupSet{i}=zeros(0,1);
            else
                dupSet{i}=(dupSet{i})*wghtSet(i); %multiply wght of daughters by wght of parent
            end
        end
        
        
        doneSet=cat(1,doneSet,nodeSet);
        doneWgt=cat(1,doneWgt,wghtSet);
        matSet=cat(1,subSet{:});

        wghtSet=cat(1,dupSet{:});
        nVG=nVG+1;
        nodeSet=mat2cell(matSet,ones(size(matSet,1),1),nVG);
    end
    doneSet=cat(1,doneSet,nodeSet);
    doneWgt=cat(1,doneWgt,wghtSet);
    
    nShapes=size(doneSet,1);
    shapeNV=zeros(nShapes,1);
    for i=1:nShapes
        shapeNV(i)=size(doneSet{i},2);
    end

    for i=2:nVG
        gSize(i).fullSet=allHoneyGraphsSize(cat(1,doneSet{shapeNV==i}),...
            cat(1,doneWgt(shapeNV==i)),G);
    end

    nGS=zeros(nVG,1);
    for i=2:nVG
        nGS(i)=size(gSize(i).fullSet,1);
    end

    exhGraphs.gSize=gSize;
    exhGraphs.nGS=nGS;
    exhGraphs.doneSet=doneSet;
    exhGraphs.doneWgt=doneWgt;
    exhGraphs.nodeSet=nodeSet;
    exhGraphs.wghtSet=wghtSet;
    exhGraphs.nVG=nVG;
    
end


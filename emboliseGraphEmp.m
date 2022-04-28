function pE=emboliseGraphEmp(G,sp,dP,nEmbInit)

    nNodes=size(G.Nodes,1);
    nEdges=size(G.Edges,1);
    
    dpv=spval(sp,rand(nEdges,1));
    
    embolism=randsample(nNodes,nEmbInit);
    untested=(1:nEdges)';

    embolised=[];
    while any(embolism)
        testEdges=[];
        testNodes=[];
        for i=1:length(embolism)
            [theseEdges,neighbours]=outedges(G,embolism(i));
            testEdges=cat(1,testEdges,theseEdges);
            testNodes=cat(1,testNodes,neighbours);
        end
        [testEdges,ia]=intersect(testEdges,untested);
        testNodes=testNodes(ia);
        [testNodes,ia]=setdiff(testNodes,union(embolised,embolism));
        testEdges=testEdges(ia);

        spread=dpv(testEdges)<dP;

        embolised=cat(1,embolised,embolism);
        embolism=testNodes(spread);
        untested=setdiff(untested,testEdges);
    end

    pE=length(embolised)/nNodes;

end
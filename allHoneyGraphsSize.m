function fullSet=allHoneyGraphsSize(nodeSets,setWgts,G)

    nSets=size(nodeSets,1);
    fullSet=struct('G',cell(nSets,1));
    for i=1:nSets
        fullSet(i).G=subgraph(G,nodeSets(i,:));
        fullSet(i).W=setWgts(i);
    end
    iso=false(nSets);
    for i=1:nSets-1
        for j=i+1:nSets
            iso(i,j)=isisomorphic(fullSet(i).G,fullSet(j).G);
            if iso(i,j)
                fullSet(i).W=fullSet(i).W+fullSet(j).W;
                fullSet(j).W=0;
            end
        end
    end
    
    fullSet=fullSet(~any(iso)');

end
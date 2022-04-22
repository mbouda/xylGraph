function margin=findMargin(G,regime)

%keyboard
%any edge that is part of 1 triangle and not 2 is a candidate outside edge
%the question becomes, what is the shortest cycle that includes its nodes
%and not the third node in the triangle?

%how many cycles 3-long is the edge part of?
    %if 2, internal
    %if 0, external
    %if 1, keep going
    
    cyc3=allcycles(G,'MaxCycleLength',3);
    cyc3=cat(1,cyc3{:});
    
    cyc4=allcycles(G,'MinCycleLength',4,'MaxCycleLength',4);
    cyc4=cat(1,cyc4{:});
    
    cyc5=allcycles(G,'MinCycleLength',5,'MaxCycleLength',5);
    cyc5=cat(1,cyc5{:});
    
    cyc6=allcycles(G,'MinCycleLength',6,'MaxCycleLength',6);
    cyc6=cat(1,cyc6{:});
    
    E=table2array(G.Edges);
    nE=size(E,1);
    marginEdge=false(nE,1);
      
    
    for i=1:nE
        
        mem3cyc=sum(ismember(cyc3,E(i,:)'),2)==2;
        n3cyc=sum(mem3cyc);
%         if n3cyc==2
%           This actually improves the algo but had to manually fix in the
%           end, so commented out for consistency
%
%             %check if third point in cyc included in other...
%             
%             cycs=cyc3(mem3cyc,:);
%             thirdA=setdiff(cycs(1,:),cycs(2,:));
%             thirdB=setdiff(cycs(2,:),cycs(1,:));
%             
%             connA=setdiff(E(any(ismember(E,thirdA),2),:),thirdA);
%             connB=setdiff(E(any(ismember(E,thirdB),2),:),thirdB);
%             
%             allConnA=all(ismember(connA,cycs(2,:)'));
%             allConnB=all(ismember(connB,cycs(1,:)'));
%             
%             if allConnA || allConnB
%                 n3cyc=1;
%             end
%         end
        if n3cyc<2
            rel3cyc=cyc3(mem3cyc,:);
            mem4cyc=sum(ismember(cyc4,E(i,:)'),2)==2;
            rel4cyc=cyc4(mem4cyc,:);
            mem5cyc=sum(ismember(cyc5,E(i,:)'),2)==2;
            rel5cyc=cyc5(mem5cyc,:);
            mem6cyc=sum(ismember(cyc6,E(i,:)'),2)==2;
            rel6cyc=cyc6(mem6cyc,:);
            
            switch regime
                case 'A'
                    notMain4=~ismember(rel4cyc,E(i,:));
                    badConn4=false(size(notMain4));
              
                    notMain5=~ismember(rel5cyc,E(i,:));
                    badConn5=false(size(notMain5));
                case 'B'
                    notMain3=~ismember(rel3cyc,E(i,:));
                    notMain4=~ismember(rel4cyc,E(i,:));

                    if size(notMain4,1)>1 && numel(rel3cyc)>0
                        v4NM=rel4cyc(notMain4);
                        v3NM=rel3cyc(notMain3);
                        isConnNM=ismember(sort(cat(2,v4NM,repmat(v3NM,[size(v4NM,1) 1])),2),E,'rows');

                        badConn4=false(size(notMain4));
                        badConn4(notMain4)=isConnNM;
                    else
                          badConn4=false(size(notMain4));
                    end

                    notMain5=~ismember(rel5cyc,E(i,:));
                    if size(notMain5,1)>1 && numel(rel3cyc)>0
                        v5NM=rel5cyc(notMain5);
                        v3NM=rel3cyc(notMain3);
                        isConnNM=ismember(sort(cat(2,v5NM,repmat(v3NM,[size(v5NM,1) 1])),2),E,'rows');

                        badConn5=false(size(notMain5));
                        badConn5(notMain5)=isConnNM;
                    else
                          badConn5=false(size(notMain5));
                    end
            end
            
            %eliminate cycles that contain entire subcycles
            bad43=sum(ismember(rel4cyc,rel3cyc) | badConn4,2)>=3;
            bad53=sum(ismember(rel5cyc,rel3cyc) | badConn5,2)>=3;
            bad63=sum(ismember(rel6cyc,rel3cyc),2)>=3;
            
            bad54=sum(ismember(rel5cyc,rel4cyc),2)>=4;
            bad64=sum(ismember(rel6cyc,rel4cyc),2)>=4;
            
            bad65=sum(ismember(rel6cyc,rel5cyc),2)>=5;
            
            bad6=bad65 | bad64 | bad63;
            bad5=bad54 | bad53;
            bad4=bad43;
            
            marginEdge(i)=n3cyc+sum(cat(1,~bad4,~bad5,~bad6))<2;
        end
    end
    
    margI=unique(E(marginEdge,:));
    margin=false(size(G.Nodes,1),1);
    margin(margI)=true;
    margin=margin | degree(G)==2;
end


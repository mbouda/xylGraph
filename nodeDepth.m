function depth=nodeDepth(G,ideal,regime)

    if ideal
        margin=degree(G)<6;
    else
        margin=findMargin(G,regime);
    end
    
    %here, impose a margin on those nodes missed out
%     margin([103 128 156 177 180 181 193 255])=true;
    
    
    D=distances(G);
    dist2bdry=D(:,margin);
    sdb=sort(dist2bdry,2,'ascend');
    if size(sdb,2)<8
        depth=geomean(sdb,2);
    else
        depth=geomean(sdb(:,1:8),2);
    end
end

% 
% dat=load('../Lycophytes/graphFiles/AdiantumPedatumGraph.mat');
% ctrInds=dat.ctrInds;
% edges=dat.edges;
% find(ctrInds(:,1)>334 & ctrInds(:,1)<336 & ctrInds(:,2)>1950 & ctrInds(:,2)<1952)
% [iE,~]=find(edges==103)
% 
% i=287
% ctrInds(edges(i,:),:)
% 
% 
% figure; plot(dat.ctrInds(:,1),dat.ctrInds(:,2),'.')
% 
% hold on
% plot(dat.ctrInds(margin,1),dat.ctrInds(margin,2),'ro')
% 
% for i=1:dat.nE
%     plot(dat.ctrInds(dat.edges(i,:),1),dat.ctrInds(dat.edges(i,:),2),'g-')
% end
% 
% for i=1:dat.nT
%     text(dat.ctrInds(i,1)+5,dat.ctrInds(i,2)-5,sprintf('%d',i))
% end

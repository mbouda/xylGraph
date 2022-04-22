function [cn,np]=recPathConcEff(orig,targ,np,cn,D)
    
    d=D(orig,targ);
    cand=find(D(:,targ)==d-1 & D(:,orig)==1);
    cn(cand,orig)=cn(cand,orig)+1;
    
    nC=size(cand,1);
    for i=1:nC
        cn(:,orig)=cn(:,orig)+cn(:,cand(i));
    end
    np(orig)=sum(np(cand));
end
function [pathConc,pc,count]=pathConcFun(G)

    D=distances(G);
    nV=size(G.Nodes,1);
    count=zeros(nV,1);
%     pc=0;
    NP=zeros(nV);
    
    for i=1:nV-1
         J=setdiff(1:nV,i)';
         [ds,ji]=sort(D(J,i),'ascend');
         jj=J(ji);

         J=jj(ds>1)';
         cn=zeros(nV);
         np=zeros(nV,1);
         np(D(:,i)==1)=1;
         for j=J
             [cn,np]=recPathConcEff(j,i,np,cn,D);
         end
         NP(:,i)=np;
         count=count+sum(cn(:,(i+1):nV),2);
%          pc=pc+sum(np((i+1):nV));
    end
    
    NP(D<=1)=0;
    pc=sum(sum(tril(NP)));
    %pathConc=sqrt(nV)*sqrt(sum(count.^2))/pc;
    pathConc=sqrt(sum(count.^2))/pc;
end
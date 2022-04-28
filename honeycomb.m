function crds=honeybcomb(nV,nL)
   
    TH=0; R=0;

    v1=[1 0];
    v2=[cos(pi/3) sin(pi/3)];
    
    for i=1:nL %18 is enough such that inceasing to 25 doesn't change order up until iVert=600
        
        dirs=[(i:-1:1)' (0:(i-1))'];
        th=zeros(i,1);
        r=zeros(i,1);
        for j=1:i
            crd=dirs(j,1)*v1+dirs(j,2)*v2;
            [th(j),r(j)]=cart2pol(crd(:,1),crd(:,2));
        end
        
        r=repmat(r,[6 1]);
        th=reshape(repmat(th,[1 6])+repmat((0:5)*pi/3,[i 1]),[6*i 1]);
        
        TH=cat(1,TH,th);
        R=cat(1,R,r);
    end
    
   % nVert=size(TH,1);
    
    [rs,ri]=sort(R);
    ths=TH(ri);
    [x,y]=pol2cart(ths,rs);
    
    crds=[x(1:nV) y(1:nV)];
    
end



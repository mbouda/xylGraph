function crds=honeyStar(nArms,width,nVMax)

    %Initialise
    TH=0; R=0;
    v1=[1 0];
    v2=[cos(pi/3) sin(pi/3)];
    v3=[sin(pi/3) cos(pi/3)];

    l0=width-1;

    %central Hex:
    for i=1:l0
        
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

    nVMin=size(TH,1);
    vIncr=nArms*width;
    nL=ceil((nVMax-nVMin)/vIncr);
    
    dirs=[(l0:-1:0)' (0:l0)'];
    L0=zeros(width,2);
    for j=1:width
        L0(j,:)=dirs(j,1)*v1+dirs(j,2)*v2;
    end
    
        
    arm=repmat((1:nL)',[1 2]).*repmat(v3,[nL 1 width])+repmat(reshape(L0',[1 2 width]),[nL 1 1]);
    arm=reshape(permute(arm,[1 3 2]),[nL*width 2]);
    [th,r]=cart2pol(arm(:,1),arm(:,2));
    
    r=repmat(r,[nArms 1]);
    
    angleLookup={[],[],[2 4],[2 3 4],[1 2 4 5],1:5};
    dTH=(2*pi/6)*angleLookup{nArms};
    
    %dTH=(1:nArms-1)*2*pi/(nArms);
    
    dTH=repmat(dTH,[nL*width 1]);
    
    th=cat(1,th,reshape(repmat(th,[1 nArms-1])+dTH,[nL*width*(nArms-1) 1]));
    
    TH=cat(1,TH,th);
    R=cat(1,R,r);
    
    %figure; polarplot(TH,R,'.')
    
    [x,y]=pol2cart(TH,R);
    
    crds=[x y];


end
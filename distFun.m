function D=distFun(x,y,nPt)

    D=sqrt((repmat(x,[1 nPt])-repmat(x',[nPt 1])).^2+...
        (repmat(y,[1 nPt])-repmat(y',[nPt 1])).^2);

end
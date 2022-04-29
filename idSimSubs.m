function match=idSimSubs(subA,subB,crds)

    nA=size(subA,1);
    nB=size(subB,1);
    nP=size(subA,2);
    
    match=false(nA,nB);
    
    ctrA=[mean(reshape(crds(subA,1),[nA nP]),2) ...
              mean(reshape(crds(subA,2),[nA nP]),2)];
    ctrB=[mean(reshape(crds(subB,1),[nB nP]),2) ...
              mean(reshape(crds(subB,2),[nB nP]),2)];
    
    THRA=cell(nA,1);
    CRDA=cell(nA,1);
    for j=1:nA
        CRDA{j}=round([crds(subA(j,:),1)-ctrA(j,1),crds(subA(j,:),2)-ctrA(j,2)],4);
        [th,r]=cart2pol(CRDA{j}(:,1),CRDA{j}(:,2));
        [rs,ri]=sort(r);
        THRA{j}=round([th(ri) rs],4);
    end

    THRB=cell(nB,1);
    CRDB=cell(nB,1);
    for j=1:nB
        CRDB{j}=round([crds(subB(j,:),1)-ctrB(j,1),crds(subB(j,:),2)-ctrB(j,2)],4);
        [th,r]=cart2pol(CRDB{j}(:,1),CRDB{j}(:,2));
        [rs,ri]=sort(r);
        THRB{j}=round([th(ri) rs],4);
    end

    for j=1:nA
        for k=1:nB
            if ~any(match(:,k))
                if all(ismember(CRDA{j},CRDB{k},'rows')) %translation
                    match(j,k)=true;
                    break  %ideally, would also not revisit subB(k,:) again...
                elseif all(THRA{j}(:,2)==THRB{k}(:,2)) %same radius

                    [thj,thji]=sort(THRA{j}(:,1));
                    [thk,thki]=sort(THRB{k}(:,1));

                    dTHJ=thj-circshift(thj,-1);
                    dTHK=thk-circshift(thk,-1);

                    dTHJ(dTHJ>pi)=dTHJ(dTHJ>pi)-2*pi; %(also -ve)
                    dTHK(dTHK>pi)=dTHK(dTHK>pi)-2*pi; %(also -ve)

                    if all(ismember(round([dTHJ THRA{j}(thji,2)],4),round([dTHK THRB{k}(thki,2)],4),'rows'))
                        match(j,k)=true;
                        break
                    else
                        thk=-thk;
                        dTHK=thk-circshift(thk,1);
                        dTHK(dTHK>pi)=dTHK(dTHK>pi)-2*pi; %(also -ve?)
                        if all(ismember(round([dTHJ THRA{j}(thji,2)],4),round([dTHK THRB{k}(thki,2)],4),'rows'))
                            match(j,k)=true;
                            break
                        end
                    end
                end
            end
        end
    end
    
end
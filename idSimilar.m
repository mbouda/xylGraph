function match=idSimilar(propSets,crds,nAdj)

        ctrs=[mean(reshape(crds(propSets,1),size(propSets)),2) ...
              mean(reshape(crds(propSets,2),size(propSets)),2)];

        THR=cell(nAdj,1);
        CRD=cell(nAdj,1);
        for j=1:nAdj
            CRD{j}=round([crds(propSets(j,:),1)-ctrs(j,1),crds(propSets(j,:),2)-ctrs(j,2)],4);
            [th,r]=cart2pol(CRD{j}(:,1),CRD{j}(:,2));
            [rs,ri]=sort(r);
            THR{j}=round([th(ri) rs],4);
        end

        match=false(nAdj);
        for j=1:nAdj-1
            for k=j+1:nAdj
                if all(ismember(CRD{j},CRD{k},'rows')) %translation
                    match(j,k)=true;
                elseif all(THR{j}(:,2)==THR{k}(:,2)) %same radius

                    [thj,thji]=sort(THR{j}(:,1));
                    [thk,thki]=sort(THR{k}(:,1));


                    dTHJ=thj-circshift(thj,-1);
                    dTHK=thk-circshift(thk,-1);

                    dTHJ(dTHJ>pi)=dTHJ(dTHJ>pi)-2*pi; %(also -ve)
                    dTHK(dTHK>pi)=dTHK(dTHK>pi)-2*pi; %(also -ve)

                    if all(ismember(round([dTHJ THR{j}(thji,2)],4),round([dTHK THR{k}(thki,2)],4),'rows'))
                        match(j,k)=true;
                    else
                        thk=-thk;
                        dTHK=thk-circshift(thk,1);
                        dTHK(dTHK>pi)=dTHK(dTHK>pi)-2*pi; %(also -ve?)
                        match(j,k)=all(ismember(round([dTHJ THR{j}(thji,2)],4),round([dTHK THR{k}(thki,2)],4),'rows'));
                    end
                end
            end
        end

end
%% This script contains the main code to drive each part of the analysis


%% Load and display graphs constructed from empirical steles


    dataDir='../Lycophytes/graphFiles/'; %directory where graphs are saved
    fileName='StenokoleosBifidusGraph.mat'; %file to be displayed
    dat=load(sprintf('%s%s',dataDir,fileName));
    
    close all
    figure(1); set(1,'Units','inches','Position',[1 1 5 5],'Color',[1 1 1])
    hA=axes('Position',[0.1 0.05 0.8 0.9]);
    imshow(dat.B)
    hold on
    plot(dat.ctrInds(:,2),dat.ctrInds(:,1),'bs','MarkerSize',2,'MarkerFaceColor','b')
    for j=1:dat.nE
        plot(dat.ctrInds(dat.edges(j,:),2),dat.ctrInds(dat.edges(j,:),1),'b-')
    end
    hold off


%% Calculate selected graph metrics
    %assumes graph loaded as dat.G, as above

    ccmp=conncomp(dat.G)';
    nComp=max(ccmp); %number separate xylem strands in image
    mDeg=mean(degree(dat.G)); %mean number neighbours
    if nComp>1    
        nVC=zeros(nComp,1);
        PCC=zeros(nComp,1);
        dpt=cell(nComp,1);
        for j=1:nComp
            nVC(j)=sum(ccmp==j);
            H=subgraph(dat.G,ccmp==j);
            PCC(j)=pathConcFun(H);
            dpt{j}=nodeDepth(H,false,'A'); %note that for very small xyelm strands (e.g. Isoetes), the count must be done manually instead as the code assumes a certain baseline size
            %The third argument for nodeDepth is 'regime' of finding the
            %margin. As the margin is not a purely topological concept,
            %finding it from just topological information (as here) is not
            %perfect. For best results, try both regimes. If needed, fix
            %manually. Actual effects found to be minimal. Depths shown in
            %Supplementary panel & used in analysis were fixed manually as
            %needed.
        end
        nVert=sum(nVC);
        pthC=sum(PCC.*nVC)/nVert; %weighted mean by number of conduits in each strand
        depth=cat(1,dpt{:}); 
        mDpt=mean(depth);
        xDpt=max(depth);
    else
        nVert=size(dat.G.Nodes,1);
        pthC=pathConcFun(dat.G);

        depth=nodeDepth(dat.G,false,'A'); %regime can be 'A' or 'B';
        mDpt=mean(depth);
        xDpt=max(depth);
    end
    nVC=nVert./nComp; %mean number of conduits in xylem strand
    
%% Synthetic graph construction I: Selected shapes

    %terete
    %strap
    %3-lobe
    %6-lobe

%% Synthetic graph construction II: Random sampling


%% Synthetic graph construction III: Constrained ellipse growth


%% Synthetic graph construction IV: Exhaustive sampling


%% Embolism spread simulation


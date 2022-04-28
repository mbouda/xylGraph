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
            %For synthetic graphs constructed with code below, marginal
            %nodes are those with fewer than 6 neighbours. Switch second
            %argument to true to make this the criterion for margin ID.
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
    
%% Load synthetic graphs (shape series in Figure 1)

    fileName=''; % options (copy into fileName string):  
                 %          hcRings / hcElls / hcTril / hcSixl
                 %files must be downloaded from OSF repository and placed
                 %on matlab path
    dat=load(fileName);


%% Synthetic graph construction I: Growth of selected shapes

    % Terete shape
     nV=120;
     nL=25;
     
     crds=honeycomb(nV,nL); %distribute nV nodes on hexagonal lattice of at most nL concentric layers
                            %the points occupy the nV innermost
                            %positions in the hexagonal lattice
     D=distFun(crds(:,1),crds(:,2),nV); 
     conn=abs(D-1)<1e-6;  %connectivity matrix: nodes separated by distance of 1 will be connected
     [ei,ej]=find(triu(conn)); 
     G=graph(ei,ej); %resulting graph
        
    % Strap shape of give width and number of ranks nr
    
        %for single file (width == 1)
        
    nV=100; %desired number of vertices

    ei=(1:nV-1)';
    ej=(2:nV)';
    G=graph(ei,ej);
    
        % for widths >1
    width=4;
    nr=10;
    
    x=1:nr;
    if width/2==ceil(width/2)
        r=width/2;
        y=(1-r):r;
    else
        r=(width-1)/2;
        y=-r:r;
    end
    [xx,yy]=meshgrid(x,y);
    xx=reshape(xx,[width*nr 1]);
    yy=reshape(yy,[width*nr 1]);
    DT=delaunayTriangulation(xx,yy);  %delaunay triangulation of points on square grid 
                                      %results in graph topologically homolgous
                                      %to one constructed using nodes on
                                      %honeycomb lattice chosen to keep
                                      %within a given width and # ranks
    E=edges(DT);
    G=graph(E(:,1),E(:,2));
    
    %Lobed shapes with central hexagon
        
    nLobes=3; %tested for 3 and 6 lobes, other values possible
    width=4;  %width of lobes & side-length of central hexagon
    nVTarg=85; %target number of nodes on graph (real value will differ to keep all lobes same length)
    
    crds=honeyStar(nLobes,width,nVTarg);
    nV=size(crds,1);  %actual number of vertices in final graph
    
    DT=delaunayTriangulation(crds(:,1),crds(:,2)); %as before, combines points
                                                   %on hex & square grids
                                                   %Delaunay triangulation
                                                   %eliminates the
                                                   %difference for this
                                                   %shape
    E=edges(DT);
    D=sqrt(sum((crds(E(:,1),:)-crds(E(:,2),:)).^2,2));
    E(D>=(sqrt(2)+1e-6),:)=[];  %but need to eliminate all edges between nonadjacent nodes
    G=graph(E(:,1),E(:,2));

    
%% Synthetic graph construction II: Random sampling
    
    nV=30; %number of vertices in desired graphs
    nG=30; %number of desired graphs
    nCores=4; %number of CPU cores available for task 
    %NOTE: should edit function below if using 1 core for small task.
    
    nPts=16.6e6;  %default value 16.6e6
    nLay=5.8e3;   %default value 5.8e3
    %Default values allow for random sampling of hex lattice on a disc of
    %radius=2140 points; 
    %reduce for memory purposes if looking for smaller graphs; need to set
    %sufficient nLay to cover desired nPts
    
    randHGS=randSampHoney(nV,nG,nCores,nPts,nLay); %returns structure containing
        %randomly constructed graphs, component arrays (coordinates, edges)
        %and metadata (rng seeds, number graphs, number vertices)

%% Synthetic graph construction III: Constrained ellipse growth

    nCores=4;
    nPts=6000; %default value 2.7e5
    nLay=60;   %default value 360   
    %These defaults allow for ra up to about 720
    
    targP88=0.62; %default value 0.62 (p88 value in sigma ASP of small terete graph)
    
    ra=16;  %choose value of major radius of ellipse in number vertices that fit
    lb=2;   %lower bound on minor radius
    ub=10;  %upper bound on minor radius

    el=driveOptConstVC(nCores,ra,lb,ub,targP88,nPts,nLay);
    %returns structure containing resulting major & minor radius, vertex coordiantes and graph


%% Embolism spread simulation

    nEmbInit=1; %number embolisms placed randomly at start of simulaiton
                %single initial embolism for standard curves in study
    asp='stdNrm';  %uses standard normal distribution as ASP
                   %set to 'emp' to use empricial distribution instead
    %iSp=0;  uncomment and set to desired species for empirical ASP

    nT=5000; %number of independent trials at each pressure
    res=64;  %resolution: number sampling points along pressure axis
    
    switch asp
        case 'stdNrm'
            X=(-4:0.01:4)';
            F=normcdf(X,0,1);
            sp=spaps(F,X,0);
            p0=-3.5;
            p1= 3.5;
            stPt=[3 0];
        case 'emp'
            aspDat=load('aspObservations.mat'); %this file needs to be placed on path
            [F,X]=ecdf(abs(aspDat.species(iSp).asp));
            X=cat(1,X(1)*0.95,X(2:end));
            sp=spaps(F,X,0);
        
            p0=spval(sp,0.015);
            p1=spval(sp,0.985);
    
            stPt=[-6 median(abs(asp.fernASP(jASP).asp))];  %start point for absolute units: [-6 & median ASP]
                 %starting point for VC fitting needs to be chosen 
                 %appropriately to the locaiton and scale of the
                 %ASP distribution. First entry should approximate the
                 %max slope and the second the median value of the
                 %expected VC
    end
    
    %the VC function runs parfor loops, assuming a parallel pool is open or
    %opening one if not. For non-parallel execution, need to edit to simple
    %for loops
    VC=findVCEmpPar(G,sp,p0,p1,nT,res,stPt,nEmbInit);
        %returns structure containing: 
        %dP values used
        %corresponding proportion conduits embolised
        %fitted VC function
        %goodness of fit statistics
        % p12, p50, & p88 values
    
        %results should always be checked for poor fit & parameters (p0,p1,stPt)
        %changed if need be.
        
        
function SetUpStaggeredMesh(dx,xmax,ymax)

global Cell

x=dx/2:dx:xmax;
y=dx/2:dx:ymax;

[X,Y]=meshgrid(x,y);

m=0;
n=0;
p=0;
q=0;
r=0;
s=0;

for i=1:size(X,1)
    for j=1:size(X,2)
        m=m+1;
        Xc(m)=X(i,j);
        Yc(m)=Y(i,j);
        ic(m)=m;
        
        if(((i~=1)&&(i~=size(X,1)))&&((j~=1)&&(j~=size(X,2))))
            n=n+1;
            ComputeCell(n)=m;        
        end 
        if(i==size(X,1))
            p=p+1;
            BdryCell{3}(p)=m;
        end
        if(i==1)
            q=q+1;
            BdryCell{4}(q)=m;
        end
        if(j==1)
            r=r+1;
            BdryCell{2}(r)=m;
        end
        if(j==size(X,2))
            s=s+1;
            BdryCell{1}(s)=m;
        end
    end
end

Pts=[Xc;Yc]';

[idx(:,1:9),~]=knnsearch(Pts,Pts(:,1:2),'K',9);
Neighbor=zeros(m,8);

for i=1:n
    index=ComputeCell(i);
    Neighbors=idx(index,2:9);
    Coordinates=[Xc(Neighbors);Yc(Neighbors)]';
    [~,numbers]=sortrows(Coordinates,[1 2]);
    Neighbor(index,1)=Neighbors(numbers(7));%East
    Neighbor(index,2)=Neighbors(numbers(2));%West
    Neighbor(index,3)=Neighbors(numbers(5));%North
    Neighbor(index,4)=Neighbors(numbers(4));%South
    Neighbor(index,5)=Neighbors(numbers(8));%NE
    Neighbor(index,6)=Neighbors(numbers(6));%SE
    Neighbor(index,7)=Neighbors(numbers(3));%NW
    Neighbor(index,8)=Neighbors(numbers(1));%SW
    

end
clear Pts;
Pts=[Xc(ComputeCell);Yc(ComputeCell)]';


for i=1:4
    BdryPtIndex=BdryCell{i}(:);
    %for j=1:size(BdryPtIndex)
        Point=BdryPtIndex;
        ProbePoint=[Xc(Point);Yc(Point)]';
        [idx,~]=knnsearch(Pts,ProbePoint,'K',1);
        Neighbor(Point(:),1)=ComputeCell(idx(:));
    %end

end

Cell.Staggered.Xc=Xc;
Cell.Staggered.Yc=Yc;
Cell.Staggered.ComputeCell=ComputeCell;
Cell.Staggered.ncellsmax=m;
Cell.Staggered.computecellsmax=n;
Cell.Staggered.Neighbor=Neighbor;
Cell.Staggered.BdryCell=BdryCell;


StaggeredCoordinates=[Xc',Yc'];


BaseCoordinates=[Cell.Xc',Cell.Yc'];


clear idx
[idx(:,1:4),~]=knnsearch(BaseCoordinates,StaggeredCoordinates,'K',4);
Cell.Base_neighbors_of_staggered=idx;

clear idx
[idx(:,1:4),~]=knnsearch(StaggeredCoordinates,BaseCoordinates,'K',4);
for i = 1:4
    Pts=Cell.BdryCell{i};
    idx(Pts,3:4)=idx(Pts,1:2);
end
clear Pts
nPtsPerDim=sqrt(size(BaseCoordinates,1));
Pts=[1,nPtsPerDim,Cell.ncellsmax-nPtsPerDim+1,Cell.ncellsmax];
idx(Pts,2)=idx(Pts,1);
idx(Pts,3:4)=idx(Pts,1:2);

Cell.Staggered_neighbors_of_base=idx;




    

    

end


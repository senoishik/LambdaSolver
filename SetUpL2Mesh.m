function SetUpL2Mesh(dx,xmax,ymax)

global Cell

x=dx/2:dx:xmax;
y=0:dx:ymax;

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

Cell.L2Mesh.Xc=Xc;
Cell.L2Mesh.Yc=Yc;
Cell.L2Mesh.ComputeCell=ComputeCell;
Cell.L2Mesh.ncellsmax=m;
Cell.L2Mesh.computecellsmax=n;
Cell.L2Mesh.Neighbor=Neighbor;
Cell.L2Mesh.BdryCell=BdryCell;


L2MeshCoordinates=[Xc',Yc'];


BaseCoordinates=[Cell.Xc',Cell.Yc'];


clear idx
[idx(:,1:2),~]=knnsearch(BaseCoordinates,L2MeshCoordinates,'K',2);
Cell.L2Mesh.Base_neighbors_of_L2Mesh=idx;

StaggeredCoordinates=[Cell.Staggered.Xc',Cell.Staggered.Yc'];

clear idx
[idx(:,1:2),~]=knnsearch(L2MeshCoordinates,StaggeredCoordinates,'K',2);


Cell.L2Mesh.L2Mesh_neighbors_of_Staggered=sort(idx,2);



end


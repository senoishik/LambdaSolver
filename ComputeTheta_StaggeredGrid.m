function ComputeTheta_StaggeredGrid()

  global Cell
  global Input
  
  InterpolateL1L2AtCellFaces()  
  GetLambdaAtStaggeredPoints()
    
  nComputeCells=Cell.Staggered.computecellsmax;
  ncellsmax=Cell.Staggered.ncellsmax;
  
  %Cell.Theta=zeros(1,ncellsmax);
  
  dx=Cell.dx;
  mu=Input.mu;
  epsilon=Input.Epsilon;
  Bm=Input.Bm;
  tbar=Input.t_bar;
 % tbar_applied=0;
  
  StiffnessMatrix=zeros(ncellsmax,ncellsmax);
  ForceVector=zeros(ncellsmax,1);
  
  ComputeCellIndex=Cell.Staggered.ComputeCell(:);
  
  %Indices
  EastIndex=Cell.L1Mesh.L1Mesh_neighbors_of_Staggered(ComputeCellIndex,2);
  WestIndex=Cell.L1Mesh.L1Mesh_neighbors_of_Staggered(ComputeCellIndex,1);
  NorthIndex=Cell.L2Mesh.L2Mesh_neighbors_of_Staggered(ComputeCellIndex,2);
  SouthIndex=Cell.L2Mesh.L2Mesh_neighbors_of_Staggered(ComputeCellIndex,1);
  
 
  
  L1=Cell.Staggered.Lambda1(ComputeCellIndex);
  L1East=Cell.L1Mesh.L1(EastIndex);
  L1West=Cell.L1Mesh.L1(WestIndex);
  
    
  
  L2=Cell.Staggered.Lambda2(ComputeCellIndex);
  L2North=Cell.L2Mesh.L2(NorthIndex);
  L2South=Cell.L2Mesh.L2(SouthIndex);
  

  Dummy = 0.0;
  dL1_dx=Central(L1,L1East,L1West,Dummy,Dummy,Dummy,Dummy,Dummy,Dummy,1,dx/2);
  dL2_dy=Central(L2,Dummy,Dummy,L2North,L2South,Dummy,Dummy,Dummy,Dummy,2,dx/2);
  
  clear EastIndex
  clear WestIndex
  clear NorthIndex
  clear SouthIndex
  
  EastIndex=Cell.Staggered.Neighbor(ComputeCellIndex,1);
  WestIndex=Cell.Staggered.Neighbor(ComputeCellIndex,2);
  NorthIndex=Cell.Staggered.Neighbor(ComputeCellIndex,3);
  SouthIndex=Cell.Staggered.Neighbor(ComputeCellIndex,4);

  
  
  
  for i=1:nComputeCells
      
      %if(CentralFlag(ComputeCellIndex(i))==1)
        StiffnessMatrix(ComputeCellIndex(i),ComputeCellIndex(i))=4.0;
        StiffnessMatrix(ComputeCellIndex(i),EastIndex(i))=-1.0;
        StiffnessMatrix(ComputeCellIndex(i),WestIndex(i))=-1.0;
        StiffnessMatrix(ComputeCellIndex(i),NorthIndex(i))=-1.0;
        StiffnessMatrix(ComputeCellIndex(i),SouthIndex(i))=-1.0; 
        ForceVector(ComputeCellIndex(i))=-1*(dL1_dx(i)+dL2_dy(i));
      %else
%         StiffnessMatrix(ComputeCellIndex(i),ComputeCellIndex(i))=60/12;
%         StiffnessMatrix(ComputeCellIndex(i),EastIndex(i))=-16/12;
%         StiffnessMatrix(ComputeCellIndex(i),WestIndex(i))=-16/12;
%         StiffnessMatrix(ComputeCellIndex(i),NorthIndex(i))=-16/12;
%         StiffnessMatrix(ComputeCellIndex(i),SouthIndex(i))=-16/12;
%         StiffnessMatrix(ComputeCellIndex(i),EastEastIndex(i))=1/12;
%         StiffnessMatrix(ComputeCellIndex(i),WestWestIndex(i))=1/12;
%         StiffnessMatrix(ComputeCellIndex(i),NorthNorthIndex(i))=1/12;
%         StiffnessMatrix(ComputeCellIndex(i),SouthSouthIndex(i))=1/12;
%         
%         L1EastEast=Cell.Staggered.Lambda1(EastEastIndex(i));
%         L1WestWest=Cell.Staggered.Lambda1(WestWestIndex(i));
%         L2NorthNorth=Cell.Staggered.Lambda2(NorthNorthIndex(i));
%         L2SouthSouth=Cell.Staggered.Lambda2(SouthSouthIndex(i));
%         
%         
%         ForceVector(ComputeCellIndex(i))=...
%          (L1EastEast-L1WestWest+L2NorthNorth-L2SouthSouth)*12/dx+...
%          (8/6)*(dL1_dx(i)+dL2_dy(i));
      %end    
  end
  StiffnessMatrix=StiffnessMatrix/(dx*dx);
  
  
  
  for i=1:4
    Pts=Cell.Staggered.BdryCell{i}(:);
    Nbrs=Cell.Staggered.Neighbor(Pts,1);
    nbdrypts=size(Pts,1);
    clear L1
    clear L2
    L1=Cell.Staggered.Lambda1(Pts);
    L2=Cell.Staggered.Lambda2(Pts);
%     if (i==1)
%         tbar_applied = tbar;
%     else
%         tbar_applied = 0;
%     end
    for j=1:nbdrypts
        StiffnessMatrix(Pts(j),Pts(j))=1.0/(dx);
        StiffnessMatrix(Pts(j),Nbrs(j))=-1.0/(dx);
        Factor1=(1-heaviside(i-2.5));
        Factor2=heaviside(i-2.5);
        ForceVector(Pts(j))=(power(-1,i+1)*tbar+power(-1,i+1)*Factor1*L1(j)+power(-1,i+1)*Factor2*L2(j));
    end
  end
  
  Corner=zeros(2,2);
  Corner(2,2)=Cell.Staggered.BdryCell{2}(1);
  Corner(1,2)=Cell.Staggered.BdryCell{1}(1);
  Corner(2,1)=Cell.Staggered.BdryCell{2}(size(Cell.Staggered.BdryCell{2},2));
  Corner(1,1)=Cell.Staggered.BdryCell{1}(size(Cell.Staggered.BdryCell{1},2));
  %Corner
  
  Xneighbor(1,1)=Corner(1,1)-1;
  Xneighbor(1,2)=Corner(1,2)-1;
  Xneighbor(2,1)=Corner(2,1)+1;
  Xneighbor(2,2)=Corner(2,2)+1;
  
  npts_per_dir=sqrt(ncellsmax);
  Yneighbor(1,1)=Corner(1,1)-npts_per_dir;
  Yneighbor(2,1)=Corner(2,1)-npts_per_dir;
  Yneighbor(2,2)=Corner(2,2)+npts_per_dir;
  Yneighbor(1,2)=Corner(1,2)+npts_per_dir;
  
  for i=1:2
      for j=1:2
          
%           if (i==1)
%                 tbar_applied = tbar;
%           else
%                 tbar_applied = 0;
%           end
          Pt=Corner(i,j);
          StiffnessMatrix(Pt,:)=0.0;
          StiffnessMatrix(Pt,Pt)=1.0/dx;
          StiffnessMatrix(Pt,Xneighbor(i,j))=-0.5/dx;
          StiffnessMatrix(Pt,Yneighbor(i,j))=-0.5/dx;
          L1=Cell.Staggered.Lambda1(Pt);
          L2=Cell.Staggered.Lambda2(Pt);
          ForceVector(Pt)=(power(-1,i+1)*0.5*tbar+0.5*((power(-1,i+1)*L1+power(-1,j+1)*L2)));
      end
  end
  
  %Prevent floating theta
  ipt=Cell.Staggered.ComputeCell(1);
  %ipt=201;
  StiffnessMatrix(ipt,:)=0;
  ForceVector(ipt)=0;
  StiffnessMatrix(ipt,ipt)=1;
   
  
  Displacement=linsolve(StiffnessMatrix,ForceVector);
  %Displacement=pcg(StiffnessMatrix,ForceVector,1e-6,30);
  Cell.Staggered.Theta=Displacement';
  
  GetThetaAtBasePoints()

end


function GetLambdaAtStaggeredPoints()

 global Cell

  BaseNeighbors=Cell.Base_neighbors_of_staggered;

  L1=Cell.Lambda1(BaseNeighbors);
  L1_Mean=mean(L1,2);
  clear L1;
  Cell.Staggered.Lambda1=L1_Mean';

  L2=Cell.Lambda2(BaseNeighbors);
  L2_Mean=mean(L2,2);
  clear L2;
  Cell.Staggered.Lambda2=L2_Mean';


end

function GetThetaAtBasePoints()

  global Cell

  StaggeredNeighbors=Cell.Staggered_neighbors_of_base;
  
  Theta=Cell.Staggered.Theta(StaggeredNeighbors);
  Theta_Mean=mean(Theta,2);
  clear Theta;
  Cell.Theta=Theta_Mean';


end

function InterpolateL1L2AtCellFaces()

  global Cell
  
  %L1 Interpolation
  BaseNeighbors=Cell.L1Mesh.Base_neighbors_of_L1Mesh;
  L1=Cell.Lambda1(BaseNeighbors);
  L1Mean=mean(L1,2);
  Cell.L1Mesh.L1=L1Mean';
  
  clear BaseNeighbors
  
  %L2 Interpolation
  BaseNeighbors=Cell.L2Mesh.Base_neighbors_of_L2Mesh;
  L2=Cell.Lambda2(BaseNeighbors);
  L2Mean=mean(L2,2);
  Cell.L2Mesh.L2=L2Mean';

end

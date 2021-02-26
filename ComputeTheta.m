function ComputeTheta()

  global Cell
  global Input
  
  
  
  nComputeCells=Cell.computecellsmax;
  ncellsmax=Cell.ncellsmax;
  
  %Cell.Theta=zeros(1,ncellsmax);
  
  dx=Cell.dx;
  mu=Input.mu;
  epsilon=Input.Epsilon;
  Bm=Input.Bm;
  tbar=Input.t_bar;
  tbar_applied=0;
  
  StiffnessMatrix=zeros(ncellsmax,ncellsmax);
  ForceVector=zeros(ncellsmax,1);
  
  ComputeCellIndex=Cell.ComputeCell(:);
  
  EastIndex=Cell.Neighbor(ComputeCellIndex,1);
  WestIndex=Cell.Neighbor(ComputeCellIndex,2);
  NorthIndex=Cell.Neighbor(ComputeCellIndex,3);
  SouthIndex=Cell.Neighbor(ComputeCellIndex,4);

  
  L1=Cell.Lambda1(ComputeCellIndex);
  L1East=Cell.Lambda1(EastIndex);
  L1West=Cell.Lambda1(WestIndex);
    
  
  L2=Cell.Lambda2(ComputeCellIndex);
  L2North=Cell.Lambda2(NorthIndex);
  L2South=Cell.Lambda2(SouthIndex);
  
  Dummy = 0.0;
  dL1_dx=Central(L1,L1East,L1West,Dummy,Dummy,Dummy,Dummy,Dummy,Dummy,1,dx);
  dL2_dy=Central(L2,Dummy,Dummy,L2North,L2South,Dummy,Dummy,Dummy,Dummy,2,dx);

  
  for i=1:nComputeCells
      StiffnessMatrix(ComputeCellIndex(i),ComputeCellIndex(i))=4.0;
      StiffnessMatrix(ComputeCellIndex(i),EastIndex(i))=-1.0;
      StiffnessMatrix(ComputeCellIndex(i),WestIndex(i))=-1.0;
      StiffnessMatrix(ComputeCellIndex(i),NorthIndex(i))=-1.0;
      StiffnessMatrix(ComputeCellIndex(i),SouthIndex(i))=-1.0;      
  end
  StiffnessMatrix=StiffnessMatrix/(dx*dx);
  ForceVector(ComputeCellIndex)=-1*(dL1_dx+dL2_dy);
  
  
  for i=1:4
    Pts=Cell.BdryCell{i}(:);
    Nbrs=Cell.Neighbor(Pts,1);
    nbdrypts=size(Pts,1);
    clear L1
    clear L2
    L1=Cell.Lambda1(Pts);
    L2=Cell.Lambda2(Pts);
    
    if (i==1)
        tbar_applied = tbar;
    else
        tbar_applied = 0;
    end
        
    for j=1:nbdrypts
        StiffnessMatrix(Pts(j),Pts(j))=1.0/(dx);
        StiffnessMatrix(Pts(j),Nbrs(j))=-1.0/(dx);
        Factor1=(1-heaviside(i-2.5));
        Factor2=heaviside(i-2.5);
        ForceVector(Pts(j))=(tbar_applied+power(-1,i+1)*Factor1*L1(j)+power(-1,i+1)*Factor2*L2(j));
    end
  end
  
  Corner=zeros(2,2);
  Corner(2,2)=Cell.BdryCell{2}(1);
  Corner(1,2)=Cell.BdryCell{1}(1);
  Corner(2,1)=Cell.BdryCell{2}(size(Cell.BdryCell{2},2));
  Corner(1,1)=Cell.BdryCell{1}(size(Cell.BdryCell{1},2));
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
          
          if (i==1)
                tbar_applied = tbar;
          else
                tbar_applied = 0;
          end
          
          Pt=Corner(i,j);
          StiffnessMatrix(Pt,:)=0.0;
          StiffnessMatrix(Pt,Pt)=1.0/dx;
          StiffnessMatrix(Pt,Xneighbor(i,j))=-0.5/dx;
          StiffnessMatrix(Pt,Yneighbor(i,j))=-0.5/dx;
          L1=Cell.Lambda1(Pt);
          L2=Cell.Lambda2(Pt);
          ForceVector(Pt)=(0.5*tbar_applied+0.5*((power(-1,i+1)*L1+power(-1,j+1)*L2)));
      end
  end
  
  %Prevent floating theta
  ipt=Cell.ComputeCell(1);
  StiffnessMatrix(ipt,:)=0;
  ForceVector(ipt)=0;
  StiffnessMatrix(ipt,ipt)=1;
   
  
  Displacement=linsolve(StiffnessMatrix,ForceVector);
  %Displacement=pcg(StiffnessMatrix,ForceVector,1e-6,30);
  Cell.Theta=Displacement';

end
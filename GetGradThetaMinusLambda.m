function [Output]=GetGradThetaMinusLambda()

global Cell

dx=Cell.dx;

Output=zeros(Cell.Staggered.ncellsmax);



ComputeCellIndex=Cell.Staggered.ComputeCell(:); 
EastIndex=Cell.Staggered.Neighbor(ComputeCellIndex,1);
WestIndex=Cell.Staggered.Neighbor(ComputeCellIndex,2);
NorthIndex=Cell.Staggered.Neighbor(ComputeCellIndex,3);
SouthIndex=Cell.Staggered.Neighbor(ComputeCellIndex,4);

Theta=Cell.Staggered.Theta(ComputeCellIndex);
ThetaEast=Cell.Staggered.Theta(EastIndex);
ThetaWest=Cell.Staggered.Theta(WestIndex);
ThetaNorth=Cell.Staggered.Theta(NorthIndex);
ThetaSouth=Cell.Staggered.Theta(SouthIndex);

  
L1=Cell.Staggered.Lambda1(ComputeCellIndex);
L2=Cell.Staggered.Lambda2(ComputeCellIndex);

Dummy=0.0;
  
dtheta_dx=Central(Theta,ThetaEast,ThetaWest,Dummy,Dummy,Dummy,Dummy,Dummy,Dummy,1,dx);
dtheta_dy=Central(Theta,Dummy,Dummy,ThetaNorth,ThetaSouth,Dummy,Dummy,Dummy,Dummy,2,dx);



Output(ComputeCellIndex)=sqrt((dtheta_dx-L1).^2+(dtheta_dy-L2).^2);


end


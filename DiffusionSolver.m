function [L1Source,L2Source]=DiffusionSolver ()


  global Cell
  global Input

  [mu1,mu2]=GetDiffusionCoefficient();


  
  nComputeCells=Cell.computecellsmax;
  Dummy=zeros(1,nComputeCells);
  dx=Cell.dx;
    
  ComputeCellIndex=Cell.ComputeCell(:);
  
  EastIndex=Cell.Neighbor(ComputeCellIndex,1);
  WestIndex=Cell.Neighbor(ComputeCellIndex,2);
  NorthIndex=Cell.Neighbor(ComputeCellIndex,3);
  SouthIndex=Cell.Neighbor(ComputeCellIndex,4);
  NEIndex=Cell.Neighbor(ComputeCellIndex,5);
  SEIndex=Cell.Neighbor(ComputeCellIndex,6);
  NWIndex=Cell.Neighbor(ComputeCellIndex,7);
  SWIndex=Cell.Neighbor(ComputeCellIndex,8);
  
  Theta=Cell.Theta(ComputeCellIndex);
  ThetaEast=Cell.Theta(EastIndex);
  ThetaWest=Cell.Theta(WestIndex);
  ThetaNorth=Cell.Theta(NorthIndex);
  ThetaSouth=Cell.Theta(SouthIndex);
  ThetaNE=Cell.Theta(NEIndex);
  ThetaSE=Cell.Theta(SEIndex);
  ThetaNW=Cell.Theta(NWIndex);
  ThetaSW=Cell.Theta(SWIndex);
  
  L1=Cell.Lambda1(ComputeCellIndex);
  L1East=Cell.Lambda1(EastIndex);
  L1West=Cell.Lambda1(WestIndex);
  L1North=Cell.Lambda1(NorthIndex);
  L1South=Cell.Lambda1(SouthIndex);
  L1NE=Cell.Lambda1(NEIndex);
  L1SE=Cell.Lambda1(SEIndex);
  L1NW=Cell.Lambda1(NWIndex);
  L1SW=Cell.Lambda1(SWIndex);
  
  
  L2=Cell.Lambda2(ComputeCellIndex);
  L2East=Cell.Lambda2(EastIndex);
  L2West=Cell.Lambda2(WestIndex);
  L2North=Cell.Lambda2(NorthIndex);
  L2South=Cell.Lambda2(SouthIndex);
  L2NE=Cell.Lambda2(NEIndex);
  L2SE=Cell.Lambda2(SEIndex);
  L2NW=Cell.Lambda2(NWIndex);
  L2SW=Cell.Lambda2(SWIndex);
                                                                                                                                                                    
  dL1_dy2=Central(L1,Dummy,Dummy,L1North,L1South,Dummy,Dummy,Dummy,Dummy,4,dx);
  dL2_dxdy=Central(L2,L2East,L2West,L2North,L2South,L2NE,L2SE,L2NW,L2SW,5,dx);
  dL1_dxdy=Central(L1,L1East,L1West,L1North,L1South,L1NE,L1SE,L1NW,L1SW,5,dx);
  dL2_dx2=Central(L2,L2East,L2West,Dummy,Dummy,Dummy,Dummy,Dummy,Dummy,3,dx);
  
   

  
  mu=Input.mu;
  epsilon=Input.Epsilon;
  tinysq=1e-32*ones(1,size(L1,2));
  absL=max(sqrt(L1.*L1+L2.*L2),tinysq);
  force_L1=0.5*(mu1/epsilon).*(2*mu/pi).*sin(2*pi*absL).*L1./absL;
  force_L2=-0.5*(mu2/epsilon).*(2*mu/pi).*sin(2*pi*absL).*L2./absL;
  %force_L1=0.0;
  %force_L2=0.0;
  
  
  %%%%%%%%%%%%%%%%%%%%%Implicit Solver%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   %Begin the Assembly
%   ncellsmax=Cell.ncellsmax;
%   StiffnessMatrix=zeros(2*ncellsmax,2*ncellsmax);
%   Force_Vector=zeros(2*ncellsmax,1);
%   
%   dsq=dx*dx;
%   delt=Cell.deltaT;
%   
%   for i=1:nComputeCells
%       
%      
%       
%       StiffnessMatrix(2*(ComputeCellIndex(i)-1)+1,2*(ComputeCellIndex(i)-1)+1)=-2.0*mu1(i)/dsq-1/delt;
%       StiffnessMatrix(2*(ComputeCellIndex(i)-1)+1,2*(NorthIndex(i)-1)+1)=1.0*mu1(i)/dsq;
%       StiffnessMatrix(2*(ComputeCellIndex(i)-1)+1,2*(SouthIndex(i)-1)+1)=1.0*mu1(i)/dsq;
%       StiffnessMatrix(2*(ComputeCellIndex(i)-1)+1,2*NEIndex(i))=-1.0*mu1(i)/dsq;
%       StiffnessMatrix(2*(ComputeCellIndex(i)-1)+1,2*SWIndex(i))=-1.0*mu1(i)/dsq;
%       StiffnessMatrix(2*(ComputeCellIndex(i)-1)+1,2*NWIndex(i))=1.0*mu1(i)/dsq;
%       StiffnessMatrix(2*(ComputeCellIndex(i)-1)+1,2*SEIndex(i))=1.0*mu1(i)/dsq;
%       
%       %Terms from the L2 Equation
%       StiffnessMatrix(2*ComputeCellIndex(i),2*ComputeCellIndex(i))=-2.0*mu1(i)/dsq-1/delt;
%       StiffnessMatrix(2*ComputeCellIndex(i),2*EastIndex(i))=1.0*mu1(i)/dsq;
%       StiffnessMatrix(2*ComputeCellIndex(i),2*WestIndex(i))=1.0*mu1(i)/dsq;
%       StiffnessMatrix(2*ComputeCellIndex(i),2*(NEIndex(i)-1)+1)=-1.0*mu1(i)/dsq;
%       StiffnessMatrix(2*ComputeCellIndex(i),2*(SWIndex(i)-1)+1)=-1.0*mu1(i)/dsq;
%       StiffnessMatrix(2*ComputeCellIndex(i),2*(SEIndex(i)-1)+1)=1.0*mu1(i)/dsq;
%       StiffnessMatrix(2*ComputeCellIndex(i),2*(NWIndex(i)-1)+1)=-1.0*mu1(i)/dsq;
%       
%       Force_Vector(2*(ComputeCellIndex(i)-1)+1,1)=-L1(i)/delt+force_L1(i);
%       Force_Vector(2*ComputeCellIndex(i),1)=-L2(i)/delt+force_L2(i);
%       
%       
%   end
%   
%   Force_Vector=-1*Force_Vector;
%   StiffnessMatrix=-1*StiffnessMatrix;
%   
%   %Apply Boundary Conditions
%   for i=1:4    
%     Pts=Cell.BdryCell{i}(:);
%     Nbrs=Cell.Neighbor(Pts,1);
%     %Cell.Lambda1(Pts)=Cell.Lambda1(Nbrs); 
%     %Cell.Lambda2(Pts)=Cell.Lambda2(Nbrs);
%     nBdryPts=length(Pts);
%     for j=1:nBdryPts
%         PtIndex=Pts(j);
%         NbrIndex=Nbrs(j);
%         StiffnessMatrix(2*(PtIndex-1)+1,2*(PtIndex-1)+1)=1.0;
%         StiffnessMatrix(2*(PtIndex-1)+1,2*(NbrIndex-1)+1)=-1.0;
%         StiffnessMatrix(2*PtIndex,2*PtIndex)=1.0;
%         StiffnessMatrix(2*PtIndex,2*NbrIndex)=-1.0;
%        
%     end
%   end
%   
%   
%   Displacement=linsolve(StiffnessMatrix,Force_Vector);
%   s=length(Displacement);
%   L1Source=Displacement(1:2:s-1);
%   L2Source=Displacement(2:2:s);
%   Cell.Lambda1=L1Source';
%   Cell.Lambda2=L2Source';
  
  
  
  %%%%%%%%%Explicit Solver%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  L1Source=mu1.*(dL1_dy2-dL2_dxdy);
  L2Source=mu2.*(dL1_dxdy-dL2_dx2);
  L1Source=L1Source-0.5*(mu1/epsilon).*(2*mu/pi).*sin(2*pi*absL).*L1./absL;
  L2Source=L2Source+0.5*(mu2/epsilon).*(2*mu/pi).*sin(2*pi*absL).*L2./absL;
  
  
  
  




end


function [Coeff1,Coeff2]=GetDiffusionCoefficient()

  global Cell
  global Input
  
  nComputeCells=Cell.computecellsmax;
  Dummy=zeros(1,nComputeCells);
  dx=Cell.dx;
  mu=Input.mu;
  epsilon=Input.Epsilon;
  Bm=Input.Bm;
  
  ComputeCellIndex=Cell.ComputeCell(:);
  
  EastIndex=Cell.Neighbor(ComputeCellIndex,1);
  WestIndex=Cell.Neighbor(ComputeCellIndex,2);
  NorthIndex=Cell.Neighbor(ComputeCellIndex,3);
  SouthIndex=Cell.Neighbor(ComputeCellIndex,4);
  NEIndex=Cell.Neighbor(ComputeCellIndex,5);
  SEIndex=Cell.Neighbor(ComputeCellIndex,6);
  NWIndex=Cell.Neighbor(ComputeCellIndex,7);
  SWIndex=Cell.Neighbor(ComputeCellIndex,8);
  
  Theta=Cell.Theta(ComputeCellIndex);
  ThetaEast=Cell.Theta(EastIndex);
  ThetaWest=Cell.Theta(WestIndex);
  ThetaNorth=Cell.Theta(NorthIndex);
  ThetaSouth=Cell.Theta(SouthIndex);
  
  
  L1=Cell.Lambda1(ComputeCellIndex);
  L1East=Cell.Lambda1(EastIndex);
  L1West=Cell.Lambda1(WestIndex);
  L1North=Cell.Lambda1(NorthIndex);
  L1South=Cell.Lambda1(SouthIndex);
  L1NE=Cell.Lambda1(NEIndex);
  L1SE=Cell.Lambda1(SEIndex);
  L1NW=Cell.Lambda1(NWIndex);
  L1SW=Cell.Lambda1(SWIndex);
  
  
  L2=Cell.Lambda2(ComputeCellIndex);
  L2East=Cell.Lambda2(EastIndex);
  L2West=Cell.Lambda2(WestIndex);
  L2North=Cell.Lambda2(NorthIndex);
  L2South=Cell.Lambda2(SouthIndex);
  L2NE=Cell.Lambda2(NEIndex);
  L2SE=Cell.Lambda2(SEIndex);
  L2NW=Cell.Lambda2(NWIndex);
  L2SW=Cell.Lambda2(SWIndex);
  
  dtheta_dx=Central(Theta,ThetaEast,ThetaWest,Dummy,Dummy,Dummy,Dummy,Dummy,Dummy,1,dx);
  dtheta_dy=Central(Theta,Dummy,Dummy,ThetaNorth,ThetaSouth,Dummy,Dummy,Dummy,Dummy,2,dx);
  dL1_dy2=Central(L1,Dummy,Dummy,L1North,L1South,Dummy,Dummy,Dummy,Dummy,4,dx);
  dL2_dxdy=Central(L2,L2East,L2West,L2North,L2South,L2NE,L2SE,L2NW,L2SW,5,dx);
  dL1_dxdy=Central(L1,L1East,L1West,L1North,L1South,L1NE,L1SE,L1NW,L1SW,5,dx);
  dL2_dx2=Central(L2,L2East,L2West,Dummy,Dummy,Dummy,Dummy,Dummy,Dummy,3,dx);
  
  dL1_dy=Central(L1,Dummy,Dummy,L1North,L1South,Dummy,Dummy,Dummy,Dummy,2,dx);
  dL2_dx=Central(L2,L2East,L2West,Dummy,Dummy,Dummy,Dummy,Dummy,Dummy,1,dx);
  
  
  curl_term=zeros(size(dL1_dy));
  for i=1:length(dL1_dy)
      curl_term(i)=dL1_dy(i)-dL2_dx(i);
%     if(abs(dL1_dy(i)-dL2_dx(i))>1e-3)
%       curl_term(i)=dL1_dy(i)-dL2_dx(i);
%     else
%       curl_term(i)=0.0;
%     end
    curl_term(i)=curl_term(i)*tanh(abs(curl_term(i))/1e-6);
  end
  %C1=-1*sign(curl_term).*(mu*(dtheta_dx-L1))/Bm;
  %C2=sign(curl_term).*(mu*(dtheta_dy-L2))/Bm;
  
   %L1Plus=L1North-L1;
   %L1Minus=L1-L1South;
   %APlus=max(sign(C1),Dummy);
   %AMinus=min(sign(C1),Dummy);
   %dL1_dy=(APlus.*L1Minus+AMinus.*L1Plus)/dx;
   
   %L2Plus=L2East-L2;
   %L2Minus=L2-L2West;
   %APlus=max(sign(C2),Dummy);
   %AMinus=min(sign(C2),Dummy);
   %dL2_dx=(APlus.*L2Minus+AMinus.*L2Plus)/dx;
   
 %   L1FluxFinal=C1.*dL1_dy-C1.*dL2_dx;
 %   L2FluxFinal=-C2.*dL1_dy+C2.*dL2_dx;
% %   
 %   Coeff1=L1FluxFinal*epsilon/Bm;
 %   Coeff2=-L2FluxFinal*epsilon/Bm;

  %curl_term=zeros(size(dL1_dy));
  %for i=1:length(dL1_dy)
  %  if(abs(dL1_dy(i)-dL2_dx(i))>1e-06)
  %    curl_term(i)=dL1_dy(i)-dL2_dx(i);
  %  else
  %    curl_term(i)=0.0;
  %  end
  %end

  Coeff1=(epsilon/Bm)*sign(curl_term).*(curl_term);
  Coeff2=-Coeff1;
end


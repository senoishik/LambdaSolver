function Allocate_and_Initialize(xmax,ymax)

global Cell

global Input

Xc=Cell.Xc;
Yc=Cell.Yc;


ncellsmax=size(Xc,2);
Lambda1=zeros(1,ncellsmax);
Lambda2=zeros(1,ncellsmax);


%Initialize Variables
c=20;
%c=0;
h=1.0;
a=2.0;

%Theta=0.0;%sin(pi*Xc/xmax).*sin(pi*Yc/ymax);
CenterX = Input.xmax/2;
CenterY = Input.ymax/2;
Radius=(Xc-CenterX).^2+(Yc-CenterY).^2;
Radius=sqrt(Radius);



Lambda2(:,((Yc>(c-h))&(Yc<(c+h))))=0.5*(1+tanh(-(Xc(:,((Yc>(c-h))&(Yc<(c+h))))-c)/a));


%Lambda1(:,((Xc>(c-h))&(Xc<(c+h))))=0.5*(1+tanh(-(Yc(:,((Xc>(c-h))&(Xc<(c+h))))-c)/a));


%Cell.Theta=Theta;
Cell.Lambda1=Lambda1;
Cell.Lambda2=Lambda2;

theta_bar=Input.theta_bar;
t_bar=Input.t_bar;

DomainBCs(theta_bar,t_bar)

ComputeTheta_StaggeredGrid()


% ComputeCellIndex=Cell.ComputeCell(:);
% EastIndex=Cell.Neighbor(ComputeCellIndex,1);
% WestIndex=Cell.Neighbor(ComputeCellIndex,2);
% NorthIndex=Cell.Neighbor(ComputeCellIndex,3);
% SouthIndex=Cell.Neighbor(ComputeCellIndex,4);
% 
% 
% Theta=Cell.Theta(ComputeCellIndex);
% ThetaEast=Cell.Theta(EastIndex);
% ThetaWest=Cell.Theta(WestIndex);
% ThetaNorth=Cell.Theta(NorthIndex);
% ThetaSouth=Cell.Theta(SouthIndex);
% L1=Cell.Lambda1(ComputeCellIndex);
% L1East=Cell.Lambda1(EastIndex);
% L1West=Cell.Lambda1(WestIndex);
%     
%   
% L2=Cell.Lambda2(ComputeCellIndex);
% L2North=Cell.Lambda2(NorthIndex);
% L2South=Cell.Lambda2(SouthIndex);
%   
% Dummy=0.0;
% dx=Cell.dx;
% 
% dL1_dx=Central(L1,L1East,L1West,Dummy,Dummy,Dummy,Dummy,Dummy,Dummy,1,dx);
% dL2_dy=Central(L2,Dummy,Dummy,L2North,L2South,Dummy,Dummy,Dummy,Dummy,2,dx);
%   
% Theta_Temp=(ThetaEast+ThetaWest+ThetaNorth+ThetaSouth)-(dL1_dx+dL2_dy)*dx*dx;
% Theta_Temp=Theta_Temp/4;
%   
% Theta=Theta_Temp;
% Cell.Theta(ComputeCellIndex)=Theta;
% 
% DomainBCs(theta_bar,t_bar)


end


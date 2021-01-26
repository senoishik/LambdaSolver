function DomainBCs(thetabar,tbar)

global Cell


dx=Cell.dx;


%Theta BCs
% for i=1:4
%     Pts=Cell.BdryCell{i}(:);
%     Cell.Theta(Pts)=thetabar;
% end

%Lambda Conditions
%for i=1:2
for i=1:4    
    Pts=Cell.BdryCell{i}(:);
    Nbrs=Cell.Neighbor(Pts,1);
    Cell.Lambda1(Pts)=Cell.Lambda1(Nbrs);%(abs(Cell.Theta(Pts)-Cell.Theta(Nbrs))/dx)-power(1,i)*tbar; 
    Cell.Lambda2(Pts)=Cell.Lambda2(Nbrs);
end

%Theta Conditions
% %East
% Pts=Cell.BdryCell{1}(:);
% Nbrs=Cell.Neighbor(Pts,1);
% Cell.Theta(Pts)=(tbar+Cell.Lambda1(Pts))*dx+Cell.Theta(Nbrs);
% 
% %West
% Pts=Cell.BdryCell{2}(:);
% Nbrs=Cell.Neighbor(Pts,1);
% Cell.Theta(Pts)=(tbar-Cell.Lambda1(Pts))*dx+Cell.Theta(Nbrs);
% 
% %North
% Pts=Cell.BdryCell{3}(:);
% Nbrs=Cell.Neighbor(Pts,1);
% Cell.Theta(Pts)=(tbar+Cell.Lambda2(Pts))*dx+Cell.Theta(Nbrs);
% 
% %South
% Pts=Cell.BdryCell{4}(:);
% Nbrs=Cell.Neighbor(Pts,1);
% Cell.Theta(Pts)=(tbar-Cell.Lambda2(Pts))*dx+Cell.Theta(Nbrs);



% for i=3:4
%     Pts=Cell.BdryCell{i}(:);
%     Nbrs=Cell.Neighbor(Pts,1);
%     Cell.Lambda2(Pts)=(abs(Cell.Theta(Pts)-Cell.Theta(Nbrs))/dx)-power(1,i)*tbar; 
%     Cell.Lambda1(Pts)=Cell.Lambda1(Nbrs);
% end


end


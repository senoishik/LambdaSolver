function AdvanceInTime()

     global Cell
     global Input
     
     thetabar=Input.theta_bar;
     tbar=Input.t_bar;
     
     time=0;
     if (Input.restart==0)
        time=0;
     elseif (Input.restart ==1)
        time= Cell.simulation_time(1);
     end
     
     T_end=Cell.T_total;
     ComputeCellIndex=Cell.ComputeCell;
     T_plot=Cell.T_plot;
     residual_rhs = 1;
     step=0;
     itimestep=0;
     while(time<=T_end && residual_rhs > 1E-2)
         [dt,WaveSpeed]=GetTimeStep();
         time=time+dt;
         itimestep=itimestep+1;
         disp(['itimestep = ',num2str(itimestep)])
         disp(['dt = ',num2str(dt)])
         disp(['residual_rhs = ',num2str(residual_rhs)])
         disp(['total time = ',num2str(time)])
         disp(['--------------------------'])
         
         
         [L1_flux,L2_flux,L1y,L2x]=GetFluxVectors();
         dL1dy=zeros(1,Cell.ncellsmax);
         dL2dx=zeros(1,Cell.ncellsmax);
         dL1dy(ComputeCellIndex)=L1y;
         dL2dx(ComputeCellIndex)=L2x;
         
         
         
         [~,L1,L2]=UpdateVariables(L1_flux,0,L2_flux,0,dt);
         Cell.Lambda1(ComputeCellIndex)=L1;
         Cell.Lambda2(ComputeCellIndex)=L2;
         DomainBCs(thetabar,tbar);
         
          [L1_Source,L2_Source]=DiffusionSolver();       
          [~,L1,L2]=UpdateVariables(0,L1_Source,0,L2_Source,dt);         
          Cell.Lambda1(ComputeCellIndex)=L1;
          Cell.Lambda2(ComputeCellIndex)=L2;
          DomainBCs(thetabar,tbar);

          %[~,~]=DiffusionSolver();
         
         ComputeTheta_StaggeredGrid()
         %Cell.sourceX(itimestep)=max(abs(L1_flux-L1_Source));
         Cell.sourceX(itimestep)=0.0;
         Cell.sourceY(itimestep)=max(abs(L2_flux-L2_Source));
         %Cell.rhs_norm(itimestep)=max(sqrt((L1_flux-L1_Source).^2 + (L2_flux-L2_Source).^2));
         Cell.rhs_norm(itimestep)=max(L2_flux-L2_Source);
         Cell.simulation_time(itimestep)=time;
         
         residual_rhs = Cell.rhs_norm(itimestep);
         
         if((rem(time,T_plot)>rem(time+dt,T_plot))||(itimestep==1))   
           step=step+1;  
           PostProcessMovieFiles(Cell.Lambda1,Cell.Lambda2,Cell.Theta,Cell.Xc,Cell.Yc,time,step,abs(dL1dy-dL2dx));
         end
          
         if(rem(itimestep,5)==0)         
             write_output_final(itimestep);
         end
         
     end
     
     write_output_final(itimestep);
     
     
end


function PostProcessMovieFiles(L1,L2,Theta,X,Y,time,step,curl)
        
     global Cell
     
     gradtheta_minus_lambda=GetGradThetaMinusLambda();
    
    
     npts=size(X,2);
     npts=sqrt(npts);
     
     
     m=0;
     for i=1:npts
         for j=1:npts
             m=m+1;
             x(i,j)=X(m);
             y(i,j)=Y(m);
             l1(i,j)=L1(m);
             l2(i,j)=L2(m);
             %gradtheta(i,j)=gradtheta_minus_lambda(m);
             theta(i,j)=Theta(m);
             Curl_Term(i,j)=curl(m);
             if(j==round(npts/2))
               line_y(i)=Y(m);
               liney_L2(i)=L2(m);
               liney_L1(i)=L1(m);
             end
             if(i==round(npts/2))
                 line_x(j)=X(m);
                 linex_L2(j)=L2(m);
                 linex_L1(j)=L1(m);
             end
         end
         
     end
     
     linex_final=line_x;
     liney_final=line_y;
     liney_L2_final=liney_L2;
     liney_L1_final=liney_L1;
     linex_L2_final=linex_L2;
     linex_L1_final=linex_L1;
     Cell.liney_L2_final(step,:)=liney_L2_final;
     Cell.liney_L1_final(step,:)=liney_L1_final;
     Cell.linex_L2_final(step,:)=linex_L2_final;
     Cell.linex_L1_final(step,:)=linex_L1_final;
     Cell.linex_final=linex_final;
     Cell.liney_final=liney_final;
     
     if(step==1)
         Cell.linex_initial=linex_final;
         Cell.liney_initial=liney_final;
         Cell.liney_L2_initial=liney_L2_final;
         Cell.liney_L1_initial=liney_L1_final;
         Cell.linex_L2_initial=linex_L2_final;
         Cell.linex_L1_initial=linex_L1_final;
     end
         linex_initial=Cell.linex_initial;
         liney_initial=Cell.liney_initial;
         liney_L2_initial=Cell.liney_L2_initial;
         liney_L1_initial=Cell.liney_L1_initial;
         linex_L2_initial=Cell.linex_L2_initial;
         linex_L1_initial=Cell.linex_L1_initial;
     %end
     
     
     
     X_staggered=Cell.Staggered.Xc; 
     Y_staggered=Cell.Staggered.Yc;
     npts_staggered=size(X_staggered,2);
     npts_staggered=sqrt(npts_staggered);
     
     
     m=0;
     for i=1:npts_staggered
         for j=1:npts_staggered
             m=m+1;
             x_staggered(i,j)=X_staggered(m);
             y_staggered(i,j)=Y_staggered(m);
             gradtheta(i,j)=gradtheta_minus_lambda(m);
         end        
     end
     
     
     
     
     linewidth=2;
     fontsize=12;
     
%      figure(1)    
%      %subplot(3,3,1)
%      contourf(x,y,l1,'LineColor','none')
%     %caxis([-0.25 0.25])
%      title(['Lambda1 at t = ',num2str(time)]);
%      xlabel('$\mathbf{X}$','Interpreter','Latex');
%      ylabel('$\mathbf{Y}$','Interpreter','Latex');
%      axis square
%      set(gca,'LineWidth',linewidth)
%      set(gca,'FontSize',fontsize)
%      colorbar
%      Cell.Lambda1_plot(step)=getframe(gcf);
     
     figure(2)
     %subplot(3,3,2)
     contourf(x,y,l2,'LineColor','none')
     %caxis([0 1.5])
     title(['Lambda2 at t = ',num2str(time)]);
     xlabel('$\mathbf{X}$','Interpreter','Latex');
     ylabel('$\mathbf{Y}$','Interpreter','Latex');
     axis square
     set(gca,'LineWidth',linewidth)
     set(gca,'FontSize',fontsize)
     colorbar
     Cell.Lambda2_plot(step)=getframe(gcf);
     
     figure(3)
     %subplot(3,3,3)
     contourf(x,y,theta,'LineColor','none')
     %caxis([-0.1 1.5])
     title(['theta at t = ',num2str(time)]);
     xlabel('$\mathbf{X}$','Interpreter','Latex');
     ylabel('$\mathbf{Y}$','Interpreter','Latex');
     axis square
     set(gca,'LineWidth',linewidth)
     set(gca,'FontSize',fontsize)
     %set(gca,'ColorScale','log')
     colorbar
     Cell.theta_plot(step)=getframe(gcf);
     
     
     figure(4)
     %subplot(3,3,3)
     contourf(x_staggered,y_staggered,(gradtheta),'LineColor','none')
     %caxis([0 0.18])
     title(['|grad(theta)-Lambda| at t = ',num2str(time)]);
     %title(['d(theta)/dy at t = ',num2str(time)]);
     xlabel('$\mathbf{X}$','Interpreter','Latex');
     ylabel('$\mathbf{Y}$','Interpreter','Latex');
     %ylim([15 25])
     axis square
     set(gca,'LineWidth',linewidth)
     set(gca,'FontSize',fontsize)
     %set(gca,'ColorScale','log')
     colorbar
     Cell.grad_theta(step)=getframe(gcf);
     
     figure(5)
     contourf(x,y,sqrt(l1.*l1+l2.*l2),'LineColor','none')
     %caxis([0 1.5])
     title(['|Lambda| at t = ',num2str(time)]);
     xlabel('$\mathbf{X}$','Interpreter','Latex');
     ylabel('$\mathbf{Y}$','Interpreter','Latex');
     axis square
     set(gca,'LineWidth',linewidth)
     set(gca,'FontSize',fontsize)
     colorbar
     Cell.abs_lambda(step)=getframe(gcf);
     
     figure(6)
     contourf(x,y,Curl_Term,'LineColor','none')
     %caxis([0 0.09])
     title(['|L1,y-L2,x| at t = ',num2str(time)]);
     xlabel('$\mathbf{X}$','Interpreter','Latex');
     ylabel('$\mathbf{Y}$','Interpreter','Latex');
     axis square
     set(gca,'LineWidth',linewidth)
     set(gca,'FontSize',fontsize)
     colorbar
     Cell.curl_term(step)=getframe(gcf);
     
     
     fontsize=8;
     figure(7)
     subplot(2,2,1)
     plot(line_y,liney_L2,'r','LineWidth',2)
     ylim([0 1.2])
     xlabel('$\mathbf{Y}$','Interpreter','Latex');
     ylabel('$\mathbf{\lambda_2}$','Interpreter','Latex');
     title(['x = 20,t = ',num2str(time)]);
     set(gca,'LineWidth',linewidth)
     set(gca,'FontSize',fontsize)
     
     
     subplot(2,2,2)
     plot(line_x,linex_L2,'r','LineWidth',2)
     ylim([0 1.2])
     xlabel('$\mathbf{X}$','Interpreter','Latex');
     ylabel('$\mathbf{\lambda_2}$','Interpreter','Latex');
     title(['y = 20, t = ',num2str(time)]);
     set(gca,'LineWidth',linewidth)
     set(gca,'FontSize',fontsize)
     
     subplot(2,2,3)
     plot(line_y,liney_L1,'r','LineWidth',2)
     ylim([-0.08 0.08])
     xlabel('$\mathbf{Y}$','Interpreter','Latex');
     ylabel('$\mathbf{\lambda_1}$','Interpreter','Latex');
     title(['x = 20,t = ',num2str(time)]);
     set(gca,'LineWidth',linewidth)
     set(gca,'FontSize',fontsize)
     
     
     subplot(2,2,4)
     plot(line_x,linex_L1,'r','LineWidth',2)
     ylim([-0.15 0.15])
     xlabel('$\mathbf{X}$','Interpreter','Latex');
     ylabel('$\mathbf{\lambda_1}$','Interpreter','Latex');
     title(['y = 20, t = ',num2str(time)]);
     set(gca,'LineWidth',linewidth)
     set(gca,'FontSize',fontsize)
     Cell.F(step)=getframe(gcf);
     
     
     figure(8)
     subplot(2,2,1)
     h(1)=plot(Cell.liney_initial,Cell.liney_L2_initial,'r','LineWidth',2,'DisplayName','Initial');
     hold on
     h(2)=plot(Cell.liney_final,Cell.liney_L2_final(end,:),'k','LineWidth',2,'DisplayName',num2str(Cell.simulation_time(end)));
     ylim([0 1.2])
     hold off
     legend(h(1:2))
     xlabel('$\mathbf{Y}$','Interpreter','Latex');
     ylabel('$\mathbf{\lambda_2}$','Interpreter','Latex');
     title(['x = 20']);
     set(gca,'LineWidth',linewidth)
     set(gca,'FontSize',fontsize)
     
     subplot(2,2,2)
     h(1)=plot(Cell.linex_initial,Cell.linex_L2_initial,'r','LineWidth',2,'DisplayName','Initial');
     hold on
     h(2)=plot(Cell.linex_final,Cell.linex_L2_final(end,:),'k','LineWidth',2,'DisplayName',num2str(Cell.simulation_time(end)));
     ylim([0 1.2])
     hold off
     legend(h(1:2))
     xlabel('$\mathbf{X}$','Interpreter','Latex');
     ylabel('$\mathbf{\lambda_2}$','Interpreter','Latex');
     title(['y = 20']);
     set(gca,'LineWidth',linewidth)
     set(gca,'FontSize',fontsize)
     
     subplot(2,2,3)
     h(1)=plot(Cell.liney_initial,Cell.liney_L1_initial,'r','LineWidth',2,'DisplayName','Initial');
     hold on
     h(2)=plot(Cell.liney_final,Cell.liney_L1_final(end,:),'k','LineWidth',2,'DisplayName',num2str(Cell.simulation_time(end)));
     ylim([-0.08 0.08])
     hold off
     legend(h(1:2))
     xlabel('$\mathbf{Y}$','Interpreter','Latex');
     ylabel('$\mathbf{\lambda_1}$','Interpreter','Latex');
     title(['x = 20']);
     set(gca,'LineWidth',linewidth)
     set(gca,'FontSize',fontsize)
     
     subplot(2,2,4)
     h(1)=plot(Cell.linex_initial,Cell.linex_L1_initial,'r','LineWidth',2,'DisplayName','Initial');
     hold on
     h(2)=plot(Cell.linex_final,Cell.linex_L1_final(end,:),'k','LineWidth',2,'DisplayName',num2str(Cell.simulation_time(end)));
     ylim([-0.15 0.15])
     hold off
     legend(h(1:2))
     xlabel('$\mathbf{X}$','Interpreter','Latex');
     ylabel('$\mathbf{\lambda_1}$','Interpreter','Latex');
     title(['y = 20']);
     set(gca,'LineWidth',linewidth)
     set(gca,'FontSize',fontsize)
     
     
        
end



function[Theta,L1,L2]=UpdateVariables(L1_flux,L1_Source,L2_flux,L2_Source,dt)
 
  global Cell
  
  
  nComputeCells=Cell.computecellsmax;
  Dummy=zeros(1,nComputeCells);
  dx=Cell.dx;
  
  
  ComputeCellIndex=Cell.ComputeCell(:);
  
  EastIndex=Cell.Neighbor(ComputeCellIndex,1);
  WestIndex=Cell.Neighbor(ComputeCellIndex,2);
  NorthIndex=Cell.Neighbor(ComputeCellIndex,3);
  SouthIndex=Cell.Neighbor(ComputeCellIndex,4);
  
  
  Theta=Cell.Theta(ComputeCellIndex);
  ThetaEast=Cell.Theta(EastIndex);
  ThetaWest=Cell.Theta(WestIndex);
  ThetaNorth=Cell.Theta(NorthIndex);
  ThetaSouth=Cell.Theta(SouthIndex);
  
  
  L1=Cell.Lambda1(ComputeCellIndex);
  L1East=Cell.Lambda1(EastIndex);
  L1West=Cell.Lambda1(WestIndex);
    
  
  L2=Cell.Lambda2(ComputeCellIndex);
  L2North=Cell.Lambda2(NorthIndex);
  L2South=Cell.Lambda2(SouthIndex);
  
  %L1=L1-dt*L1_flux+dt*L1_Source;
  L1=0;
  L2=L2-dt*L2_flux+dt*L2_Source;
  
  

end



function [L1FluxFinal,L2FluxFinal,dL1_dy,dL2_dx]=GetFluxVectors()

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
  
  
  dL1_dy=Central(L1,Dummy,Dummy,L1North,L1South,Dummy,Dummy,Dummy,Dummy,2,dx);
  dL2_dx=Central(L2,L2East,L2West,Dummy,Dummy,Dummy,Dummy,Dummy,Dummy,1,dx);
  
  curl_term=zeros(size(dL1_dy));
  for i=1:length(dL1_dy)
    curl_term(i)=dL1_dy(i)-dL2_dx(i);
    %if(abs(dL1_dy(i)-dL2_dx(i))>1e-3)
    %  curl_term(i)=dL1_dy(i)-dL2_dx(i);
    %else
    %  curl_term(i)=0.0;
    %end
    curl_term(i)=curl_term(i)*tanh(abs(curl_term(i))/1e-06);
  end
  max(abs(curl_term))
  
  
  C1=-sign(curl_term).*(mu*(dtheta_dx-L1))/Bm;
  C2=sign(curl_term).*(mu*(dtheta_dy-L2))/Bm;
  
%   L1Plus=L1North-L1;
%   L1Minus=L1-L1South;
%   APlus=max(-C1,Dummy);
%   AMinus=min(-C1,Dummy);
%   L1Flux=(APlus.*L1Minus+AMinus.*L1Plus)*dx;
%   L2Plus=L2East-L2;
%   L2Minus=L2-L2West;
%   APlus=max(C1,Dummy);
%   AMinus=min(C1,Dummy);
%   L2Flux=(APlus.*L2Minus+AMinus.*L2Plus)*dx;
%   L1FluxFinal=L1Flux+L2Flux;
%   
%   L1Plus=L1North-L1;
%   L1Minus=L1-L1South;
%   APlus=max(-C2,Dummy);
%   AMinus=min(-C2,Dummy);
%   L1Flux=(APlus.*L1Minus+AMinus.*L1Plus)*dx;
%   L2Plus=L2East-L2;
%   L2Minus=L2-L2West;
%   APlus=max(C2,Dummy);
%   AMinus=min(C2,Dummy);
%   L2Flux=(APlus.*L2Minus+AMinus.*L2Plus)*dx;
%   L2FluxFinal=L1Flux+L2Flux;

   L1Plus=L1North-L1;
   L1Minus=L1-L1South;
   L2Plus=L2East-L2;
   L2Minus=L2-L2West;
   
   APlus=max(C1,Dummy);
   AMinus=min(C1,Dummy);
   dL1_dy=(APlus.*L1Minus+AMinus.*L1Plus)/dx;
   APlus=max(C1,Dummy);
   AMinus=min(C1,Dummy);
   dL2_dx=(APlus.*L2Minus+AMinus.*L2Plus)/dx;
   L1FluxFinal=dL1_dy-dL2_dx;
   
   
   APlus=max(C2,Dummy);
   AMinus=min(C2,Dummy);
   dL1_dy=(APlus.*L1Minus+AMinus.*L1Plus)/dx;
   APlus=max(C2,Dummy);
   AMinus=min(C2,Dummy);
   dL2_dx=(APlus.*L2Minus+AMinus.*L2Plus)/dx;
   L2FluxFinal=-dL1_dy+dL2_dx;
   
   
   
   
   
   
 end

function [dt,WaveSpeed]=GetTimeStep()

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
  
  dtheta_dx=Central(Theta,ThetaEast,ThetaWest,Dummy,Dummy,Dummy,Dummy,Dummy,Dummy,1,dx);
  dtheta_dy=Central(Theta,Dummy,Dummy,ThetaNorth,ThetaSouth,Dummy,Dummy,Dummy,Dummy,2,dx);
  dL1_dy2=Central(L1,Dummy,Dummy,L1North,L1South,Dummy,Dummy,Dummy,Dummy,4,dx);
  dL2_dxdy=Central(L2,L2East,L2West,L2North,L2South,L2NE,L2SE,L2NW,L2SW,5,dx);
  dL1_dxdy=Central(L1,L1East,L1West,L1North,L1South,L1NE,L1SE,L1NW,L1SW,5,dx);
  dL2_dx2=Central(L2,L2East,L2West,Dummy,Dummy,Dummy,Dummy,Dummy,Dummy,3,dx);
  
  dL1_dy=Central(L1,Dummy,Dummy,L1North,L1South,Dummy,Dummy,Dummy,Dummy,2,dx);
  dL2_dx=Central(L2,L2East,L2West,Dummy,Dummy,Dummy,Dummy,Dummy,Dummy,1,dx);
  
%   Non_Convex(:,1)=(2*mu/pi)*sin(2*pi*abs(L1+L2)).*L1./sqrt(L1.*L1+L2.*L2)
%   Non_Convex(:,2)=(2*mu/pi)*sin(2*pi*abs(L1+L2)).*L2./sqrt(L1.*L1+L2.*L2)
  
  WaveSpeed(1,:)=sign(dL1_dy-dL2_dx).*(-mu*(dtheta_dx-L1))/Bm;%-epsilon*(dL1_dy2-dL2_dxdy))/Bm;
  WaveSpeed(2,:)=sign(dL1_dy-dL2_dx).*(-mu*(dtheta_dy-L2))/Bm;%+epsilon*(dL1_dxdy-dL2_dx2))/Bm;
  
  MaxSpeed=max(max(abs(WaveSpeed)));
  dt=0.1*dx/MaxSpeed;
  Cell.deltaT=dt;
  

end




close all; clear all; clc

global Cell

global Input

global Solution

ReadInputVariables()

SetupBaseMesh(Cell.dx,Input.xmax,Input.ymax)
SetUpStaggeredMesh(Cell.dx,Input.xmax,Input.ymax)
SetUpL1Mesh(Cell.dx,Input.xmax,Input.ymax)
SetUpL2Mesh(Cell.dx,Input.xmax,Input.ymax)


Allocate_and_Initialize(Input.xmax,Input.ymax)

AdvanceInTime()



% writerObj=VideoWriter('L1.avi');
% F=Cell.Lambda1_plot;
% writerObj.FrameRate = 1;
% % open the video writer
% open(writerObj);
% % write the frames to the video
% for i=1:length(F)
%     % convert the image to a frame
%     frame = F(i) ;    
%     writeVideo(writerObj, frame);
% end
% % close the writer object
% close(writerObj);

writerObj=VideoWriter('L2.avi');
F=Cell.Lambda2_plot;
writerObj.FrameRate = 1;
% open the video writer
open(writerObj);
% write the frames to the video
for i=1:length(F)
    % convert the image to a frame
    frame = F(i) ;    
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);

writerObj=VideoWriter('Theta.avi');
F=Cell.theta_plot;
writerObj.FrameRate = 1;
% open the video writer
open(writerObj);
% write the frames to the video
for i=1:length(F)
    % convert the image to a frame
    frame = F(i) ;    
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);


writerObj=VideoWriter('grad_theta_minus_lambda.avi');
F=Cell.grad_theta;
writerObj.FrameRate = 1;
% open the video writer
open(writerObj);
% write the frames to the video
for i=1:length(F)
    % convert the image to a frame
    frame = F(i) ;    
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);


writerObj=VideoWriter('abs_lambda.avi');
F=Cell.abs_lambda;
writerObj.FrameRate = 1;
% open the video writer
open(writerObj);
% write the frames to the video
for i=1:length(F)
    % convert the image to a frame
    frame = F(i) ;    
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);


writerObj=VideoWriter('curl_term.avi');
F=Cell.curl_term;
writerObj.FrameRate = 1;
% open the video writer
open(writerObj);
% write the frames to the video
for i=1:length(F)
    % convert the image to a frame
    frame = F(i) ;    
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);

writerObj=VideoWriter('LinePlots.avi');
F=Cell.F;
writerObj.FrameRate = 1;
% open the video writer
open(writerObj);
% write the frames to the video
for i=1:length(F)
    % convert the image to a frame
    frame = F(i) ;    
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);


linewidth=3;
fontsize=12;
     
figure(9)
plot(Cell.simulation_time,Cell.sourceX,'r','LineWidth',linewidth)
xlabel('$\mathbf{t}$','Interpreter','Latex');
ylabel('$\mathbf{RHS_1}$','Interpreter','Latex');
set(gca,'LineWidth',linewidth)
set(gca,'FontSize',fontsize)

figure(10)
plot(Cell.simulation_time,Cell.sourceY,'r','LineWidth',linewidth)
xlabel('$\mathbf{t}$','Interpreter','Latex');
ylabel('$\mathbf{RHS_2}$','Interpreter','Latex');
set(gca,'LineWidth',linewidth)
set(gca,'FontSize',fontsize)

figure(11)
plot(Cell.simulation_time,Cell.rhs_norm,'r','LineWidth',linewidth)
xlabel('$\mathbf{t}$','Interpreter','Latex');
ylabel('$\mathbf{norm of residual (rhs)}$','Interpreter','Latex');
set(gca,'LineWidth',linewidth)
set(gca,'FontSize',fontsize)

filename = 'equilibrium_feb22.mat';
save(filename)




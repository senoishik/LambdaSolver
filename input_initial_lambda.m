function input_initial_lambda()

global Cell

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

Lambda2(:,((Yc>(c-h))&(Yc<(c+h))))=0.5*(1+tanh(-(Xc(:,((Yc>(c-h))&(Yc<(c+h))))-c)/a));
%Lambda1(:,((Xc>(c-h))&(Xc<(c+h))))=0.5*(1+tanh(-(Yc(:,((Xc>(c-h))&(Xc<(c+h))))-c)/a));

%Cell.Theta=Theta;

Cell.Lambda1=Lambda1;
Cell.Lambda2=Lambda2;

end


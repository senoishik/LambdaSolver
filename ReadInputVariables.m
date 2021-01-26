function ReadInputVariables()


global Input
global Cell

A=load('Input.txt');

Cell.dx=A(1);
Input.xmax=A(2);
Input.ymax=A(3);
Input.Epsilon=A(4);
Input.mu=A(5);
Input.Bm=A(6);

Input.theta_bar=A(7);
Input.t_bar=A(8);

Cell.T_total=A(9);
Cell.T_plot=A(10);



end


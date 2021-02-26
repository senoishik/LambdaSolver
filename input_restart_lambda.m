function input_restart_lambda()

global Cell

L1 = load('lambda1.txt');
L2 = load('lambda2.txt');

Cell.simulation_time(1)=L1(1);
Cell.Lambda1 = L1(2:end)';
Cell.Lambda2 = L2(2:end)';

end


function write_output_final(step)

global Cell

str1= ['outputs/lambda1-',num2str(step),'.txt'];
str = join(str1);
fileID1= fopen(str,'w');
fprintf(fileID1,'%12.16f\r\n',Cell.simulation_time(1,end));
fprintf(fileID1,'%12.16f\r\n',Cell.Lambda1');
fclose(fileID1);

str1= ['outputs/lambda2-',num2str(step),'.txt'];
str = join(str1);
fileID2= fopen(str,'w');
fprintf(fileID2,'%12.16f\r\n',Cell.simulation_time(1,end));
fprintf(fileID2,'%12.16f\r\n',Cell.Lambda2');
fclose(fileID2);


end


function [Derivative]=Central(X,East,West,North,South,NE,SE,NW,SW,type,dx)

   switch type
       
       case 1 %d/dx
           Derivative=(East-West)/(2*dx);
       case 2 %d/dy
           Derivative=(North-South)/(2*dx);    
       case 3 %d2/dx2
           Derivative=(East+West-2*X)/(dx*dx);
       case 4 %d2/dy2
            Derivative=(North+South-2*X)/(dx*dx);
       case 5 %d2/dxdy    
            Derivative=(SW+NE-SE-NW)/(4*dx*dx);
       
       
   end
   
   

end
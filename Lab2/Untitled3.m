i = 1;

   for Qi =  0.1:0.3:300
       Ri=0.0002*Qi
       flag = Lab2_CD_const(Ri,Qi);
       disp(Ri);
       disp(Qi);
       if(flag == 1)
          RQ(:,i) = [Ri;Qi];
          i = i+1;
          
       end
       
   end
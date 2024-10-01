function y = func_test(a)
    t = 0;
   parfor i = 1:30
     t(i) = a;
     if(t(i) == 20)
         error('a too large');
     end
   end 
   y = t;
end


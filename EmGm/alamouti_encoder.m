function [tx1,tx2] = alamouti_encoder(data)
   data = (1/(2)^0.5)*data;
   for ii = 1:length(data)
       if(mod(ii,2)==1)
           tx1(ii,1) = data(ii);
           tx2(ii,1) = data(ii+1);
       else
           tx1(ii,1) = -conj(data(ii));
           tx2(ii,1) = conj(data(ii-1));
       end
   end
end




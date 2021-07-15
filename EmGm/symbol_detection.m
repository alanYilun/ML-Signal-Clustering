function [s1,S] = symbol_detection(y,Ref)
 %% symbol detection
    x = size(y,2);
    for jj = 1:length(y(:,1))
        for mm = 1:16
        if y(jj,x) == mm
           s1(jj,1:2) = Ref(mm,3:4);
           s1(jj,3) = y(jj,x);
        end
        end
    end
    
    % Alamouti Decoder
    for k = 1:length(s1)
        S(2*k-1,1) = s1(k,1);
        S(2*k,1) = s1(k,2);
        S(2*k-1,2) = s1(k,3);
        S(2*k,2) = s1(k,3);
    end
    
end
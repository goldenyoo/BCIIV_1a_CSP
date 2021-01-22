function [check, result] = myClassifier(fp)
    load('C:\Users\유승재\Desktop\true_labels\feature.mat');
    
      
    check = (0.5*(fp - Mr)'*pinv(Qr)*(fp - Mr) + 0.5*log(det(Qr))) - (0.5*(fp - Ml)'*pinv(Ql)*(fp - Ml) + 0.5*log(det(Ql))); 
 
    if check > 0
        result = -1;  % l
    else
        result = 1;   % r
    end
end

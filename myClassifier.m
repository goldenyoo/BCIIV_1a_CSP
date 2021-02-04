function [check, result] = myClassifier(fp,M1,M2,Q1,Q2)
    
      
    check = (0.5*(fp - M1)'*pinv(Q1)*(fp - M1) + 0.5*log(det(Q1))) - (0.5*(fp - M2)'*pinv(Q2)*(fp - M2) + 0.5*log(det(Q2))); 
 
    if check > 0
        result = -1;  % 2
    else
        result = 1;   % 1
    end
end

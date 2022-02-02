function tis=tisentropic(Tt1,pt2,pt1,k);
% Tt1= total temperature at initial state
% pt2 = total pressure at final state
% pt1 = total pressure at initial state
% k = cp/cv (average Temperature)
    exp=(k-1)/k;
    tis=Tt1*((pt2/pt1)^exp);
end
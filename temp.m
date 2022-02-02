%calculate temperature from total temperature
function T=temp(Tt,c,cp)
    T=Tt-c^2/(2*cp);
end
function lValue = CalcLVector(xvalue,tvalue)

    lValue = zeros(1,length(tvalue));
    for i=1:length(tvalue(:,1))
        lValue(i) = (CalculateL(xvalue(i,1)) - CalculateLStar());
    end

end
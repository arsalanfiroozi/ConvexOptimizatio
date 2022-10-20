function y = phi(S,a,b,g)
    y = 0;
    for i=1:length(S)
        y = y + a + b*S(i) + g*S(i)^2;
    end
end


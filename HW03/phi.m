function z = phi(x,y)
    z = 0;
    eps = 0.1;
    for i=1:length(x)
        z = z + (x(i) + y(i))^2 + eps * (x(i)^2+y(i)^2);
    end
end


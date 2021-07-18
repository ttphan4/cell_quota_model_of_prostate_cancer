function x = u(t)
global change
num = sum(change <= t);
if mod(num,2) == 0
    x = 1;
    
else
    x = ~1;
end
end



function output = binomial(n,k)

x = 1;
y = n;
for i = 1:k-1
    
    y = y*(n-i);
    x = x*(i+1);  

end

output = y/x;
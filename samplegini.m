%function to calculate the Gini from the sample.

function g = samplegini(x)

n = size(x,1);

y = sort(x);

tot = sum(y);
sumfirst = 0;
sumsecond = 0;
for i = 1:n
    f(i,1) = i/n;
end

for i = 1:n
    eta(i,1) = sum(y(1:i))/tot;
end
%calculate the first term
for i = 1:(n-1)
    
    first = eta(i+1,1)*f(i,1);
    sumfirst = sumfirst + first;
    
end

for i = 1:(n-1)
    
    second = eta(i,1)*f(i+1,1);
    sumsecond =sumsecond +second;
    
end

g = sumfirst-sumsecond;

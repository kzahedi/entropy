function[ output ] = movavg(x,n)

for t=n+1:length(x)-n; output(t) = mean(x(t-n:t+n)); end

end

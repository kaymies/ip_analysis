figure; loglog(f,x1x,'o')
one = polyfit(log(f),log(x1x),1)
hold on
loglog(f,x2xstiff,'x')
loglog(f,x4xstiff,'+')

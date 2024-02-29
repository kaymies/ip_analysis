f = @(x,k1,k2,k3,k4) sqrt(abs((k1 - k4)^2 + 4.*(2.*x+k3)*(k2-2.*x)'));
xs = -300:10:300;
plot(xs,f(xs,Gains(1,1,2),Gains(1,2,2),Gains(2,1,2),Gains(2,2,2)),'r')
xlabel('x'); ylabel('f(x)')
fun= zeros(length(xs),1);
for i = 1:length(xs)
    fun(i)= f(xs(i),Gains(1,1,2),Gains(1,2,2),Gains(2,1,2),Gains(2,2,2))
end
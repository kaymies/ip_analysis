ol = [1.3127 0.4570 -1.3127 -0.4569]';
b5 = [-1.3136
   -1.3118
   -0.4570
   -0.4570]';
b4 = [   -1.3135
   -1.3119
   -0.4570
   -0.4570]';
b3 = [  -1.3134
   -1.3120
   -0.4570
   -0.4570]';
b2 = [-1.3133
   -1.3121
   -0.4570
   -0.4570]';
b02 = [ -1.3132
   -1.3123
   -0.4570
   -0.4570]';

figure;
hold on
plot(real(ol),imag(ol),'x','MarkerSize',20,'LineWidth',3,'MarkerEdgeColor','k')

plot(real(b5),imag(b5),'+','MarkerSize',15,'LineWidth',3,'MarkerEdgeColor','r')
plot(real(b4),imag(b4),'+','MarkerSize',10,'LineWidth',3,'MarkerEdgeColor','g')
plot(real(b3),imag(b3),'+','MarkerSize',10,'LineWidth',3,'MarkerEdgeColor','m')
plot(real(b2),imag(b2),'+','MarkerSize',10,'LineWidth',3,'MarkerEdgeColor','c')
plot(real(b02),imag(b02),'+','MarkerSize',10,'LineWidth',3,'MarkerEdgeColor','b')
legend('open loop', 'b = 5', 'b = 4', 'b = 3', 'b = 2', 'b = 0.2')
xlabel('Real [Hz]')
ylabel('Imaginary [Hz]')
set(gca,'FontSize',20)
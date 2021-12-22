%n=1.25+7.25i;%aluminium
%n=1.5+8i;
n=0.07+4.2i;%silver
n=0.2+3.44i;%silver
n=0.09+3i;

o_i=(0.01:0.01:pi/2);
o_t=asin( 1/n*sin(o_i) );
R=-cos(o_i-o_t)./cos(o_i+o_t);
P=abs(R);
delta=angle(R);

figure(10);
plot(o_i/pi*180,delta,'.-');
%xlim([60 80])
figure(11)
plot(o_i/pi*180,P,'.-');
%xlim([60 80])
function g=nlGain(x)

% if x<=0
%      a=5658;
%      b=132.85;
% %g=132*ones(length(x(:,1)),length(x(1,:)));
%  g=((-b+sqrt(b^2+4*a*x))/(2*a)./x).^(-1);
% 
% 
% else
% %     a=5658;
% %     b=92;
%     g=100*ones(length(x(:,1)),length(x(1,:)));
% end
% 
a=5658;
b=132.85;
g=((-b+sqrt(b^2+4*a*x))/(2*a)./x).^(-1);

%g=1*ones(length(x(:,1)),length(x(1,:)));


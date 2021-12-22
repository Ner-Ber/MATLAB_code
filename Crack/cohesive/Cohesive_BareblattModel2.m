function [hOfx, rho] = Cohesive_BareblattModel2(X,L,d,w0,p0,LL)
    if nargin<6
        LL = 1;
    end
    %% from Landau Lifshitz p.129
    D = 5e8;
    if LL
        
        p = @(x) 0.*x + p0;
        %         g = @(x)  my_Heaviside(x-L+d).*my_Heaviside(-(x-L)).*((x-L+d).^2).*w0./(d^2);
        %         g = @(x) my_Heaviside(x-L+d).*my_Heaviside(-(x-L)).*((w0*15*pi)/16).*(x - (L-d)).^2;
        %     g = @(x)  heaviside(x-L+d).*(x-L+d).*w0./d;
        %     g = @(x)  heaviside(x-L+d).*50*w0;
        g = @(x)  my_Heaviside(x).*my_Heaviside(-(x-d)).*w0;
%         g = @(x)  my_Heaviside(x).*w0.*exp(-abs(x/d));
        omega = @(x)  (p(x)-g(x))./D;
        I = @(x,z) (omega(x))./(sqrt(L^2-z.^2).*(z+x-L));
        %--- solveIntegrals
        IofX = nan(size(X));
        for i = 1:length(X)
            xi = X(i);
            Itmp = @(z) I(xi,z);
            IofX = integral(Itmp,0,L);
        end
        rho = 2.*pi^(-2).*sqrt(X.*(2*L-X)).*IofX;
        
        
        
        %         hOfx = trapz(X,rho)-cumtrapz(X,rho);
        hOfx = cumtrapz(X,rho);
        %         hOfx = nan(size(X));
        %         for i = 1:(length(X)-1)
        %             hOfx(i) = trapz(X(i:end),rho(i:end));
        %         end
        
        
        %% from Broberg  p.95
    else
        p = @(x) 0.*x + 0; %p0;
        %         g = @(x)  my_Heaviside(x).*my_Heaviside(-(x-1)).*((x-1).^2).*w0;
%                 g = @(x)  my_Heaviside(x).*my_Heaviside(-(x-1)).*w0;
%                 g = @(x) my_Heaviside(x).*my_Heaviside(-(x-1)).*((15*pi)/16).*(-x + 1).^2;
%                 g = @(x) my_Heaviside(x).*my_Heaviside(-(x-1)).*(w0).*(-x + 1).^2;
        %         g = @(x) my_Heaviside(x).*my_Heaviside(-(x-1)).*x.^-0.5;
        %         g = @(x) real(w0*sqrt(0.5^2-(x-0.5).^2));
        %     g = @(x)  heaviside(x-L+d).*(x-L+d).*w0./d;
        %     g = @(x)  heaviside(x-L+d).*50*w0;
        %         g = @(x)  my_Heaviside(x).*my_Heaviside(-(x-1)).*(1-x).*w0 + my_Heaviside(-x).*my_Heaviside((x+1)).*(1+x).*w0;
                g = @(x)  my_Heaviside(x).*w0.*exp(-abs(x/d));
%         g = @(x) w0*my_Heaviside(x).*real(sqrt(1^2-x.^2));
        omega = @(x)  p(x)-g(x);
        I = @(e,gamma) omega(gamma).*log(abs((sqrt(e)+sqrt(gamma))./(sqrt(e)-sqrt(gamma))));   % integrant in Broberg 3.6.12
        %         figure; plot(X,g(X)); title(['g(x), \omega_0=',num2str(w0)]);
        %--- solveIntegrals
        IofX = nan(size(X));
        for i = 1:length(X)
            ei = X(i);
            Itmp = @(gamma) I(ei,gamma);
            IofX(i) = integral(Itmp,0,1);
        end
        
        hOfx = sqrt(X)-(1./2*pi).*IofX;
        
        rho = [];
    end
end
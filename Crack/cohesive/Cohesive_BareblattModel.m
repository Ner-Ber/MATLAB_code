function hOfx = Cohesive_BareblattModel(X,L,d,w0,p0)
    LL = 1;
    %% from Landau Lifshitz p.129
    D = 5e8;
    if LL
        
        p = @(x) 0.*x + p0;
        %         g = @(x)  my_Heaviside(x-L+d).*my_Heaviside(-(x-L)).*((x-L+d).^2).*w0./(d^2);
%         g = @(x) my_Heaviside(x-L+d).*my_Heaviside(-(x-L)).*((w0*15*pi)/16).*(x - (L-d)).^2;
        %     g = @(x)  heaviside(x-L+d).*(x-L+d).*w0./d;
        %     g = @(x)  heaviside(x-L+d).*50*w0;
        g = @(x)  my_Heaviside(x).*w0.*exp(-abs(x/d));
        omega = @(x)  (p(x)-g(x))./D;
        I = @(x,z) (omega(x))./(sqrt(L^2-z.^2).*(z-x));
        %--- solveIntegrals
        epsilon = 1e-13;
        IofX = nan(size(X));
        for i = 1:length(X)
            xi = X(i);
            Itmp = @(z) I(xi,z);
            if xi<-L || xi>L
                I1 = integral(Itmp,-L,L);
                I2=0;
            elseif xi==-L || xi==L
                I1 = integral(Itmp,-L+epsilon,L-epsilon);
                I2=0;
            else
                I1 = integral(Itmp,-L,xi-epsilon);
                I2 = integral(Itmp,xi+epsilon,L);
            end
            IofX(i) = I1 + I2;
%             IofX = integral(Itmp,-L,L);
        end
        rho = -pi^(-2).*sqrt(L^2-X.^2).*IofX;
        
        
        
        hOfx = trapz(X,rho)-cumtrapz(X,rho);
        %         hOfx = nan(size(X));
        %         for i = 1:(length(X)-1)
        %             hOfx(i) = trapz(X(i:end),rho(i:end));
        %         end
        
        
        %% from Broberg  p.95
    else
        p = @(x) 0.*x + 0; %p0;
        %         g = @(x)  my_Heaviside(x).*my_Heaviside(-(x-1)).*((x-1).^2).*w0;
        %         g = @(x)  my_Heaviside(x).*my_Heaviside(-(x-1)).*w0;
        %         g = @(x) my_Heaviside(x).*my_Heaviside(-(x-1)).*((15*pi)/16).*(-x + 1).^2;
        %         g = @(x) my_Heaviside(x).*my_Heaviside(-(x-1)).*x.^-0.5;
        %         g = @(x) real(w0*sqrt(0.5^2-(x-0.5).^2));
        %     g = @(x)  heaviside(x-L+d).*(x-L+d).*w0./d;
        %     g = @(x)  heaviside(x-L+d).*50*w0;
%         g = @(x)  my_Heaviside(x).*my_Heaviside(-(x-1)).*(1-x).*w0 + my_Heaviside(-x).*my_Heaviside((x+1)).*(1+x).*w0;
        g = @(x)  my_Heaviside(x).*w0.*exp(-abs(x/d));
        %         omega = @(x)  (p(x)-g(x))./D;
        omega = @(x)  g(x);
        I = @(e,gamma) omega(gamma).*log(abs((sqrt(e)+sqrt(gamma))./(sqrt(e)-sqrt(gamma))));   % integrant in Broberg 3.6.12
%         figure; plot(X,g(X)); title(['g(x), \omega_0=',num2str(w0)]);
        %--- solveIntegrals
        %         epsilon = 1e-9;
        IofX = nan(size(X));
        for i = 1:length(X)
            ei = X(i);
            Itmp = @(gamma) I(ei,gamma);
            %         if sqrt(ei)<0 || sqrt(ei)>1
            %             I1 = integral(Itmp,0,1);
            %             I2=0;
            %         elseif sqrt(ei)==0 || sqrt(ei)==1
            %             I1 = integral(Itmp,0+epsilon,1-epsilon);
            %             I2=0;
            %         else
            %             I1 = integral(Itmp,0,sqrt(ei)-epsilon);
            %             I2 = integral(Itmp,sqrt(ei)+epsilon,1);
            %         end
            %         IofX(i) = I1 + I2;
            IofX(i) = integral(Itmp,0,1);
        end
        
        hOfx = sqrt(X)-(1./2*pi).*IofX;
    end
end
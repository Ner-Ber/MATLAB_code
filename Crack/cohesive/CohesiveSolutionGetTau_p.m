function tau_p = CohesiveSolutionGetTau_p(v,Gamma,L,varargin)
    % tau_p = CohesiveSolutionGetTau_p(v,Gamma,L,cohesiveFun,PlaneFlag,x_range,dx)
    %
    % CohesiveSolutionGetTau_p will calculate value of tau_p for a cohesive
    % zone model given Gamma and L (cohesive length).
    %
    % INPUTS:
    %v - velocity in units of [Cr]
    %Gamma - fracture energy [J/m^2]
    %L - cohesive length
    %
    % OPTIONAL:
    % cohesiveFun - function describing the cutoff of the cohesive zone:
    %               options: 'exp' (default), 'linear', 'Ohnaka89', 'slipWeak'
    %               ***NOTICE*** all 'tau' hadnles should obey to the
    %               consitions tau(0)=1 and tau(-1)~0. otherwise this will work
    %               but with  a discontinuity.
    % PlaneFlag - define as 'PlaneStrain' (default) or 'PlaneStress'
    % x_range - spacial vector definning solution range (defaults: 80E-3:-dx:-80E-3 [m])
    %           Also can be a cosume x vector, that has to obey
    %           length(x_range)>2.
    % dx - spacial vector x spacing. (defaults: 5e-6 [m])
    %
    % all formulas taken from Samudrala@Samudrala,Huang&Rosakis JGR2002
    
    %% set defaults and parameters
    warning('off');
    dx_def = 5e-6;   %[m]
    x_range_def = [80E-3,-80E-3]; %[m] distance from the crack tip ,theta=0;
    [cohesiveFun,PlaneFlag,x_range,dx] = setDefaults4function(varargin,'exp','PlaneStrain',x_range_def,dx_def);
    
    if length(x_range)>2
        x = x_range;
    else
        x = max(x_range):-abs(dx):min(x_range); %[m] distance from the crack tip ,theta=0;
    end
    
    [Cd, Cs, Cr, nu , ~, ~, mu, GammaDefault, PlaneStrain, tau_p_def]=CrackSolutionMaterialProperties(PlaneFlag);
    
    k=Cs/Cd;%Broberg p.330
    v=v*Cr;
    
    %------- General LEFM functions
    alpha_d=(1-(v/Cd).^2).^0.5;
    alpha_s=(1-(v/Cs).^2).^0.5;
    D=4*alpha_d.*alpha_s-(1+alpha_s.^2).^2; % the Rayleight function
    
    if ~strcmpi(PlaneFlag,'PlaneStrain')
        A=2*(1-k^2)*alpha_s.*v.^2./D/Cs^2; % this is plain stress %Broberg p.334,336, equation 6.2.26
    else
        A=alpha_s.*v.^2./D/(1-nu)/Cs^2; %Only for plain strain
    end
    
    
    %% Create \tau profile
    
    if strcmpi(cohesiveFun,'linear') %-- Linear cohesive zone
        tau=@(x) 1+x;
        
    elseif strcmpi(cohesiveFun,'slipWeak') %-- Linear slip weakening
        tau=@(x)1-(-x).^1.4;
        
    elseif strcmpi(cohesiveFun,'exp') %-- exponential cohesive
        tau=@(x) exp((x./0.1));
        
    elseif strcmpi(cohesiveFun,'Ohnaka89') %---- Cohesive model from Ohnaka&Yamashita JGR 89
        %     tau=@(x) (1+log(1-10*(x./10))).*exp((x./10));
        tau=@(x) (1+log(1-100*x)).*exp(10*x);
        
    elseif strcmpi(cohesiveFun,'Const') %---- step function
        tau=@(x) my_Heaviside(-x).*my_Heaviside(x+1);
        
    elseif strcmpi(cohesiveFun,'MySigmoid') %---- step function
        tau = @(x) (1 - x.^2).*(0.5*tanh(3*x + 1) + 0.5) + (x + 1).^2.*(0.5*tanh(-3*x - 1) + 0.5);
        
    elseif isa(cohesiveFun,'function_handle')
        tau= cohesiveFun;
        
    elseif loadCohesive %-------
        % tmp=load('C:\Users\owner\Documents\MATLAB\tmp\c4.mat');
        % p=tmp.p;
        % % tmp=load('C:\Users\owner\Documents\MATLAB\tmp\c2.mat');
        % % p2=tmp.p;
        % ILowerBoundery=-1;
        % tau=@(x) fnval(p,x);
    end
    
    
    %% Calc tau_p
    %Note that unlike eq.13 L is not necessary ILowerBoundery.
    
    K=(Gamma*mu*4*(1-k^2)./A).^0.5;%Broberg p.334,336 Eq. 6.2.42
    tauDecay=@(x) tau(x).*(-x).^-0.5;
    %     tau_p=K/(2/pi)^0.5/(integral(tauDecay,-1,0)*L^0.5);%Eq.16c after zeta/L->x
    tau_p=K*(sqrt(L*2/pi).*integral(tauDecay,-1,0)).^-1;    %Eq.16c after zeta/L->x
   
    warning('on');
end
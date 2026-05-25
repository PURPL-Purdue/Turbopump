function [mdot,Re] = laby_leakage(c, s, w, R, N, Pi, Pe, rho, mu)


%c = 5.08e-5; % input('Enter radial clearance in m : ');
%s = .0010414; %input('Enter tooth pitch in m: ');
%w = .000254; %input('Enter tooth width in m: ');
%R = .0085; %input('Enter shaft radius in m: ');
%N = 5; %input('Enter number of teeth : ');
%Pi = 6.895e+5; %input('Enter seal inlet pressure in Pa: ');
%Pe = 544738; %input('Enter seal exit pressure in Pa: ');
%rho= 800; % input('Enter density in kg/m^3 : ');
%mu = 1.92e-3; %input('Enter dynamic viscosity in PaS : ');
A = 2*3.142*R*c;
p = ones(1,N+1); p(1) = Pi; p(N+1) = Pe;
mdot = 0.1*A*sqrt((Pi-Pe)*rho);
error=1;
for i=1:10000
 Re = mdot/(3.142*2*R*mu);
 gamma =((1-6.5*(c/s))- 8.638*(c/s) *(w/s))*(Re+((1-6.5*(c/s))- 8.638*(c/s) *(w/s))^(-1/(2.454*(c/s)+2.258*(c/s)*(w/s)^1.673)))^(2.454*(c/s)+2.258*(c/s)*(w/s)^1.673);
 cd1 =(1.097-0.0029*w/c)*(1+(44.86*w/c)/Re)^(-0.2157)/sqrt(2);
 tdc =cd1*0.925*gamma^.861;
 p(2) = p(1)-(mdot/A)^2/(2*cd1^2*rho);

 for j = 2:N-1
 p(j+1) =p(j)-(mdot/A)^2/(2*tdc^2*rho);
 end
 mdot1 = A* tdc*sqrt(2*(rho*(p(N)-p(N+1))));
 error = abs((mdot1-mdot)/mdot);
 mdot=(mdot1-mdot)*0.1+mdot;
 if(error<0.0001)
 break;
 end
end
Q = mdot/rho;

%fprintf('leakage rate in kg/s is : %f\n', mdot);
%fprintf("leakage rate in m^3/s is %f\n", Q);
%disp('The pressure distribution is (in Pa) : ');
%p 

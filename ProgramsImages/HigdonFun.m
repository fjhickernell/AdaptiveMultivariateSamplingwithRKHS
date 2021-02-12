function y = HigdonFun(x)
% HIGDON (2002) FUNCTION
% f(s) = sin(2Pi s/10)+0.2sin(2Pi s/2.5)
% s \in [0, 10]
% Let s = 10x, x \in [0,1] 

y = sin(2*pi*x)+0.2*sin(8*pi*x);

%Some other candidate functions
%y = x.*sin(10*x); 
%y = -sin(20*pi*x)./(4*x+1)+(2*x-0.5).^4; %GRAMACY & LEE (2012) FUNCTION
%y = 1 - exp(-1./(2*x)); %CURRIN ET AL. (1988) SURVIVAL FUNCTION
%y = log(1 + 5*x); %CURRIN ET AL. (1988) SURVIVAL FUNCTION
end

% Small investigation of the Plasma sheath ode in one dimension. Imaginary
% parts seem to crop around about twice the debye length. Some pretty
% interesting behavior on altering the inital conditions as well.
% Prgrmr: Hayden Orleth-Diener
% Date: 1.30.2021
tRange = [0 1];
X0 = [0;
      1];
Xo = [0;
     0.0001];
  
  
% Plasma sheath with variable M
subplot(2,2,1)
[eSol, XSol] = ode15s(@(e,X)plasma_sheath(e,X,1),tRange,Xo);
x = XSol(:,1);
plot(eSol,x)


Mn = linspace(1e-2,1e3,10);
Mnstr = double2str(Mn,'Mn');

hold on

for n = 1:length(Mn)
    [eSol, XSol] = ode15s(@(e,X)plasma_sheath(e,X,Mn(n)),tRange,Xo);
    x = XSol(:,1);
    plot(eSol,x)
end
hold off
legend(Mnstr,'Location','best')
title('Plasma Sheath with variable $$ m = u_{0}/\sqrt{KT_{e}/M} $$',...
        'Interpreter','Latex')
xlabel('$$ \xi = x/\lambda_{D} $$','Interpreter','Latex')
ylabel('$$ \chi = -e\phi/KT_{e} $$','Interpreter','Latex')

% Focuses M plasma sheath
subplot(2,2,2)
[eSol, XSol] = ode15s(@(e,X)plasma_sheath(e,X,1),tRange,Xo);
x = XSol(:,1);
plot(eSol,x)


Mn = linspace(1,1e3,10);
Mnstr = double2str(Mn,'Mn');
axis([0.7497864 0.7497875 -0.5938724 -0.5938709])
hold on


for n = 1:length(Mn)
    [eSol, XSol] = ode15s(@(e,X)plasma_sheath(e,X,Mn(n)),tRange,Xo);
    x = XSol(:,1);
    plot(eSol,x)
end
hold off
legend(Mnstr(2:end))
title('Plot 1 focused at range $$ M_{n}  = [112,1000] $$','Interpreter',...
        'Latex')
xlabel('$$ \xi = x/\lambda_{D} $$','Interpreter','Latex')
ylabel('$$ \chi = -e\phi/KT_{e} $$','Interpreter','Latex')

%Plasma sheath with variable dxde condition
subplot(2,2,3)
[eSol2, XSol2] = ode15s(@(e,X)plasma_sheath(e,X,1),tRange,X0);
x = XSol2(:,1);
plot(eSol2,x)
hold on
X2n = logspace(0,-6,7);
for j = 1:length(X2n)
    X0prime = [1;
               X2n(j)];
    [eSol2, XSol2] = ode15s(@(e,X)plasma_sheath(e,X,100),tRange,X0prime);
    x = XSol2(:,1);
    plot(eSol2,x)
end
X2nstr = double2str(X2n,'d\chi_{0}/d\xi');
legend(X2nstr,'Location','northeastoutside')
title('Plasma Sheath with variable initial condition $$ d\chi_{0}/d\xi $$',...
        'Interpreter','Latex')
xlabel('$$ \xi = x/\lambda_{D} $$','Interpreter','Latex')
ylabel('$$ \chi = -e\phi/KT_{e} $$','Interpreter','Latex')
hold off

% Plasma sheath with variable x0 condition
subplot(2,2,4)
tRange1 = [0 10];
[eSol, XSol] = ode15s(@(e,X)plasma_sheath(e,X,1),tRange,X0);
x = XSol(:,1);
plot(eSol,x)
hold on
X1n = linspace(1e-2,1e2,6);
X1nstr = double2str(X1n,'\chi_{0}');
for k = 1:length(X1n)
    X1prime = [X1n(k);
               1];
    [eSol3, XSol3] = ode15s(@(e,X)plasma_sheath(e,X,1),tRange1,X1prime);
    x = XSol3(:,1);
    plot(eSol3,x)
end
legend(X1nstr,'location','best')
title('Plasma Sheath with variable inital condition $$ \chi_{0} $$',...
        'Interpreter','Latex')
xlabel('$$ \xi = x/\lambda_{D} $$','Interpreter','Latex')
ylabel('$$ \chi = -e\phi/KT_{e} $$','Interpreter','Latex')
hold off


%function to convert a number array to string cell array
function str = double2str(double_vector,str_name)
    str = cell(1,length(double_vector));
    for n = 1:length(double_vector)
        str{1,n} = [str_name, ' = ', num2str(double_vector(n))];
    end    
end

% function for the Plasma sheath ODE
function dXde = plasma_sheath(e,X,M)
    x = X(1);
    Ve = X(2);
    dxde = Ve;
    dVede = varyM(M);
    function dVede = varyM(m)
        dVede = -1/sqrt((1 + (2*x/m^2))) - exp(-x);
    end
    dXde = [dxde;
           dVede];
end
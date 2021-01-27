%Model for a charged particle Larmor motion for perpendicular E and B
%fields. The parameters q/M*E and Larmor frequency are varied
%logarythmically to show variation in flight path.
%programmer: Hayden Orleth-Diener
%date: 1/26/21
tRange = [0 2];
Y0 = [0;
      0;
      0;
      1;
      0;
      0];

%Initial solution with constant parameters.
[tSol YSol] = ode45(@(t,Y)larmorMotion(t,Y,0.1,1),tRange,Y0);
x = YSol(:,1);
y = YSol(:,2);
z = YSol(:,3);
plot3(x,y,z)

c_n = logspace(-2,-4,6);%varied parameter q/M*E
b_n = logspace(2,3,6);  %variable parameter omegaC
cnbn = [c_n b_n];
cnbn_str = {};
for j = 1:length(cnbn)
    cnbn_str = [cnbn_str, num2str(cnbn(j))];
end
hold on

%Group of graphs over c_n
for n = 1:length(c_n)
    [tSol YSol] = ode45(@(t,Y)larmorMotion(t,Y,c_n(n),100),tRange,Y0);
    x = YSol(:,1);
    y = YSol(:,2);
    z = YSol(:,3);
    plot3(x,y,z)
    
end

%Group of graphs over b_n
for i = 1:length(b_n)
    [tSol YSol] = ode45(@(t,Y)larmorMotion(t,Y,c_n(1),b_n(n)),tRange,Y0);
    x = YSol(:,1);
    y = YSol(:,2);
    z = YSol(:,3);
    plot3(x,y,z)
end
legend(cnbn_str)
xlabel('x')
ylabel('y')
zlabel('z')
title('Larmor motion for variable Q/M*E and Larmor frequency')
hold off

%Larmor motion coupled ode
function dYdt = larmorMotion(t,Y,c0,b0)
    x = Y(1);
    y = Y(2);
    z = Y(3);
    v_x = Y(4);
    v_y = Y(5);
    v_z = Y(6);
    dxdt = v_x;
    dydt = v_y;
    dzdt = v_z;
    [dv_xdt,dv_ydt] = vary_cb(c0,b0);
    function [dv_xdt,dv_ydt] = vary_cb(c,b)
        dv_xdt = c + v_y;
        dv_ydt = b*-v_x;
    end
    dv_zdt = 1;
    dYdt = [dxdt;
            dydt;
            dzdt;
            dv_xdt;
            dv_ydt;
            dv_zdt];
end
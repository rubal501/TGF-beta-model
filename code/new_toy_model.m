%%%% new model %%%%

tgfb = 1;
n1 = 2;
k1 = 1;
gtr = 0.2;
k2 = 1;
dse = 0.2; 
as7 = 0.5;
ds7 = 0.2;
dtr = 0.2;
atr = 0.01;
ase = 0.01;
tspan = [0 100];
y0 = [0.1 0 0];

[t, y] = ode45(@(t,x)tgf(t,x,tgfb,n1,k1,gtr,k2,dse,as7,ds7,dtr,atr,ase),tspan,y0);

figure (1)
plot(t,y(:,1))
hold on
plot(t,y(:,2))
plot(t,y(:,3))
hold off

%%
syms TGFr smadE smadR
dx1 = atr + ((tgfb*TGFr^n1)/(k1^n1+TGFr^n1))-gtr*smadR*TGFr-dtr*TGFr;
dx2 = ase + (1/(1+smadR))*((TGFr*smadE)/(k2+smadE))-dse*smadE;
dx3 = as7*smadE-ds7*smadR;

[soltr,solsme,solsmr] = vpasolve(dx1==0,dx2==0,dx3==0, [TGFr smadE smadR]);

%%
syms TGFr smadE smadR tgfb n1 k1 gtr k2 dse as7 ds7 dtr atr ase

dx1 = atr + ((tgfb*TGFr^n1)/(k1^n1+TGFr^n1))-gtr*smadR*TGFr-dtr*TGFr;
dx2 = ase + (1/(1+smadR))*((TGFr*smadE)/(k2+smadE))-dse*smadE;
dx3 = as7*smadE-ds7*smadR;

vars = [TGFr,smadE,smadR];
cx1   = solve(dx1 == 0, smadR);
cx2 = solve(dx2 == 0, smadR);    
cx3  = solve(dx3 == 0, smadR);

cx1 =(atr - TGFr*dtr + (TGFr^n1*tgfb)/(TGFr^n1 + k1^n1))/(TGFr*gtr);
cx2 = - (TGFr*smadE)/((k2 + smadE)*(ase - dse*smadE)) - 1;
cx3 = (as7*smadE)/ds7;

params = struct('tgfb',1,'n1',2,'k1',1,'gtr',0.2,'k2' ,1,'dse',0.2, 'as7', ...
    0.5,'ds7',0.2,'dtr',0.2,'atr',0.01,'ase',0.01);
f1_num = subs(cx1, params);
f2_num = subs(cx2, params);
f3_num = subs(cx3, params);

f1_func = matlabFunction(f1_num, 'Vars', [TGFr, smadE, smadR]);
f2_func = matlabFunction(f2_num, 'Vars', [TGFr, smadE, smadR]);
f3_func = matlabFunction(f3_num, 'Vars', [TGFr, smadE, smadR]);

range = linspace(0, 50, 200);  
[TR,SE,SR] = meshgrid(range, range, range);

F1_vals = f1_func(TR,SE,SR);
F2_vals = f2_func(TR,SE,SR);
F3_vals = f3_func(TR,SE,SR);

figure;
p1 = patch(isosurface(TR,SE,SR, F1_vals, 0));
set(p1, 'FaceColor', [0.85 0.1 0.1], 'EdgeColor', 'none', 'FaceAlpha', 0.6); hold on;

p2 = patch(isosurface(TR,SE,SR, F2_vals, 0));
set(p2, 'FaceColor', [0.5 0.5 1], 'EdgeColor', 'none', 'FaceAlpha', 0.6);

p3 = patch(isosurface(TR,SE,SR, F3_vals, 0));
set(p3, 'FaceColor', [0 0.4 0], 'EdgeColor', 'none', 'FaceAlpha', 0.6);

xlabel('TGFr'); ylabel('smadE'); zlabel('smadR');
legend([p1 p2 p3], {'f1 = 0 (TGFr)', 'f2 = 0 (smadE)', 'f3 = 0 (smadR)'});
title('Nulclinas');
view(3); axis tight; grid on;
camlight headlight; lighting gouraud;
rotate3d on;

%%
function dx = tgf(t,x,tgfb,n1,k1,gtr,k2,dse,as7,ds7,dtr,atr,ase)

TGFr = x(1);
smadE = x(2);
smadR = x(3);

dx1 = atr + ((tgfb*TGFr^n1)/(k1^n1+TGFr^n1))-gtr*smadR*TGFr-dtr*TGFr;
dx2 = ase + (1/(1+smadR))*((TGFr*smadE)/(k2+smadE))-dse*smadE;
dx3 = as7*smadE-ds7*smadR;

dx = [dx1;dx2;dx3];

end
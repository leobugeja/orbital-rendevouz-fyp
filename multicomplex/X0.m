function X0 = intial_state(x,y,z,xdot,ydot,zdot)

x = multicomplex(inputconverter(x,1,h));
x = num2cell(x.zn);

y = multicomplex(inputconverter(y,2,h));
y = num2cell(y.zn);

z = multicomplex(inputconverter(z,3,h));
z = num2cell(z.zn);

xdot = multicomplex(inputconverter(xdot,4,h));
xdot = num2cell(xdot.zn);

ydot = multicomplex(inputconverter(ydot,5,h));
ydot = num2cell(ydot.zn);

zdot = multicomplex(inputconverter(zdot,6,h));
zdot = num2cell(zdot.zn);

X0 = {x,y,z,xdot,ydot,zdot}';
end
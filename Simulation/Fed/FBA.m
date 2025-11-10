
vg_all=[10 9.5 9 8.7 8.5 8 7.5 7 6.5 6 5.5 5.2 5 4.5 4.2 4];
vo_all=[15 14 13 12 11 10 9 8 7 6 5 4 3 2 1 0];

%E.coli core index
eth=48  ;
obj=25  ;
glu=52  ;
o2=60  ;
ac=44;
ATP=16;

%Reading the model
% model=csvread('e_coli_core.csv');
model=load('e_coli_core.mat');
model=model.e_coli_core;

A=[];
b=[];
f=-model.c;
c_z=zeros(size(model.c));
f_max = c_z;
f_max(ac)=-1;
f_min = c_z;
f_min(ac)=1;
Aeq=model.S;
beq=model.b;
lb=model.lb;
ub=model.ub;
v0=zeros(size(lb,1),1);
lb(ac)=0.0;
% lb(ATP)=10.0;


for i =1:16

for j= 1:16
%constraints
lb(obj)=0;
ub(obj)=1000;
lb(glu)=-vg_all(i);
lb(o2)=-vo_all(j);

%FBA
%Linear
[v] = linprog(f,A,b,Aeq,beq,lb,ub);
va_FBA(i,j)=v(ac);
mi=v(obj);

%FVA
lb(obj)=mi*0.95;
ub(obj)=mi*1.05;
[v_max] = linprog(f_max,A,b,Aeq,beq,lb,ub);
va_max(i,j)=v_max(ac);
[v_min] = linprog(f_min,A,b,Aeq,beq,lb,ub);
va_min(i,j)=v_min(ac);

end
end

surf(vg_all,vo_all,va_max-va_min)
xlabel('glucose uptake')
ylabel('oxygen uptake')
zlabel('Max Acetate flux - Min acetate flux')
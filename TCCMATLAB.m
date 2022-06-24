%parâmetros do sistema
m = 0.111;
Raio = 0.05;
g = -9.8;
L = 1.3;
d = 0.03;
J = 9.99e-6;
fc = 0.02;
x0=0;
y0=0;

%Representação por Estados
A=[0 1;0 -fc/(J/Raio^2+m)];
B=[0;-m*g*d/(L*(J/Raio^2+m))];
C=[1 0];
D=[0];
SistemaMA = ss(A,B,C,D);

%Polos do sistema
syms s
polos=solve(s^2+(0.690107*2.8981*s)+2.8981^2==0,s);
p1= -1.5+1.5i;
p2= -1.5-1.5i;
p3= -7.5;

%Controlabilidade e Observabilidade
Cm = [B A*B];
Om = [C;C*A];
postoC= rank(Cm);
postoO= rank(Om);

%Encontrar Ganhos por alocação de polos
Atil = [0 1 0;0 -fc/(J/Raio^2+m) 0;-1 0 0];
Btil = [B;0];
K = place(Atil,Btil,[p1 p2 p3]);
Kt = [K(1) K(2)];
Ki = -K(3);

%Encontrar Ganhos por LQR
Q=[4 0;0 0.1];
R=[0.01];
Klqr=lqr(A,B,Q,R);
s = size(A,1);
Z = [zeros([1,s]) 1];
N = inv([A,B;C,D])*Z';
Nx = N(1:s); Nu = N(1+s);
Kbar=Nu + Klqr*Nx;

%Equações realimentadas
Arob=[A-B*Kt B*Ki;-C 0];
Brob=[0;0;1];
Crob=[1 0 0];
Sistemarobusto=ss(Arob,Brob,Crob,D);
SistemaLQR=ss(A-B*Klqr,B*Kbar,C,D);

%Gráfico à Ganho por polos
t = 0:0.01:10;
u = 1*ones(size(t));
lsim(Sistemarobusto,u,t);

%Gráfico à Ganho LQR
r = ones(size(t));
lsim(SistemaLQR,r,t);


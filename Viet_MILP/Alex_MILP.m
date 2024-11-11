% Canonical Form
% min(x) f(T).x  f(T) la ma tran chuyen vi cua F
% s.t
% x(intcon) are integers
% A.x <= b
% Aeq.x = beq
% lb <= x <= ub
% whereas x, f, b, beq, lb, ub are vectors. x = [x1;x2;...xn]
%         A, Aeq are matrices

% Problem:
% min(x) x1 + 5x2 + 6x3 - 2x4
% s.t
% x3, x4 are integers
% x1 + 4x2 >= -5 <=> -x1 - 4x2 <= 5
% x2 - x4 <= 5
% x1 + 3x3 - 10 x4 <= 10
% x1 + 2x2 - 3x3 + 4x4 <= 2
% -20 <= x1, x2 <= 20
% 0 <= x3, x4 <= 1

clear, pack, clc
f = [1;5;6;-2];
intcon = [3;4];
b = [5;5;10;2];
A = [-1 -4 0 0;
    0 1 0 -1;
    1 0 3 -10;
    1 2 -3 4];

Aeq = [];
beq = [];
lb = [-20; -20; 0; 0];
ub = [20; 20; 1; 1];
x0 = []; %% 

options = optimoptions('intlinprog','Display','iter');
[sol,fval,exitflag,output] = intlinprog(f,intcon,A,b,Aeq,beq,lb,ub,x0,options)

%% simpler example
clear, pack, clc
f=[1; 5; 6; -2];
intcon=[3;4]; % specify optimization variables that are integers
b=[5;5;10;2];
A=[-1 -4 0 0; 
    0 1 0 -1; 
    1 0 3 -10; 
    1 2 -3 4];
Aeq=[];
beq=[];
lb=-20*ones(4,1);
ub=20*ones(4,1);
x0=[];
 
options = optimoptions('intlinprog','Display','iter');
 
[sol,fval,exitflag,output] = intlinprog(f,intcon,A,b,Aeq,beq,lb,ub,x0,options)




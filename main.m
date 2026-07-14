clc;
clear;
close all;

% ------------------------

X = [
0;
4;
6;
2;
3
];

Y = [
0;
0;
0;
3;
5
];

% ------------------------

Conn=[
1 2;
2 3;
1 4;
2 4;
3 4;
4 5;
3 5
];

% ------------------------

N = length(X); % Number of Nodes
M = size(Conn, 1); % Number of Members

% ------------------------

E=[
200e9;
200e9;
70e9;
200e9;
200e9;
70e9;
200e9
];


A=[
0.001;
0.001;
0.002;
0.001;
0.001;
0.002;
0.001
];


% ------------------------
F_ext=[
0 0;
0 0;
0 0;
0 0;
0 -10000
];

% ------------------------
% Boundary Conditions
% 1 = Fixed
% 0 = Free

BC=[
1 1;
1 1;
0 0;
0 0;
0 0
];


% ------------------------
% Global Stiffness Matrix

K_total = zeros(2*N,2*N);


for m = 1:M

    i = Conn(m,1);
    j = Conn(m,2);

    L = sqrt((X(j)-X(i))^2 + (Y(j)-Y(i))^2);

    C = (X(j)-X(i))/L;
    S = (Y(j)-Y(i))/L;


    k_element = E(m)*A(m)/L*[
        C^2    C*S   -C^2   -C*S;
        C*S    S^2   -C*S   -S^2;
       -C^2   -C*S    C^2    C*S;
       -C*S   -S^2    C*S    S^2];


    index = [2*i-1 2*i 2*j-1 2*j];


    K_total(index,index)=K_total(index,index)+k_element;

end


% ------------------------
% Find Free DOFs

free_dofs = [];

for node = 1:N

    if BC(node,1)==0

        dof_x = 2*node-1;
        free_dofs(end+1)=dof_x;

    end


    if BC(node,2)==0

        dof_y = 2*node;
        free_dofs(end+1)=dof_y;

    end

end


disp("Free DOFs:")
disp(free_dofs)


% ------------------------
% Reduced System

K_reduced = K_total(free_dofs,free_dofs);



% Convert External Force to Vector

F = zeros(2*N,1);


for node=1:N

    F(2*node-1)=F_ext(node,1);
    F(2*node)=F_ext(node,2);

end


F_reduced = F(free_dofs);



% ------------------------
% Check Stability

if rank(K_reduced)<size(K_reduced,1)

    error("Structure is unstable. Singular stiffness matrix.")

end



% ------------------------
% Solve

Delta_reduced = K_reduced\F_reduced;


Delta=zeros(2*N,1);


for i=1:length(free_dofs)

    Delta(free_dofs(i))=Delta_reduced(i);

end



% ------------------------
% Print Displacement

fprintf("\nDisplacement Results\n")
fprintf("---------------------\n")

for node=1:N

    fprintf("Node %d : Ux = %.8e , Uy = %.8e\n",...
        node,...
        Delta(2*node-1),...
        Delta(2*node));

end



% ------------------------
% Reaction Forces

R = K_total*Delta-F;

fprintf("\nReaction Forces\n")
fprintf("---------------------\n")


for node = 1:N
    R(abs(R)<1e-8)=0;
    fprintf("Node %d : Rx = %.3f N , Ry = %.3f N\n",...
        node,...
        R(2*node-1),...
        R(2*node));
end

% ------------------------
% Member Forces

fprintf("\nMember Forces\n")
fprintf("--------------------\n")

for m = 1:M

    i = Conn(m, 1);
    j = Conn(m, 2);

    L = sqrt((X(j)-X(i))^2 + (Y(j)-Y(i))^2);

    C = (X(j)-X(i))/L;
    S = (Y(j)-Y(i))/L;

    d = [
        Delta(2*i-1);
        Delta(2*i);
        Delta(2*j-1);
        Delta(2*j)];

    N_force = E(m)*A(m)/L * [-C -S C S]*d;
    Stress = N_force / A(m);

    if abs(N_force)<1e-6

        type="Zero Force";

    elseif N_force > 0

        type="Tension";

    else

        type="Compression";

    end
    fprintf("Member %d : N = %.3f N , Stress = %.3e Pa (%s)\n",...
        m,...
        N_force,...
        Stress,...
        type)
end
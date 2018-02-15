clear all;
clc;
format long
%Inputs: Interior points
m=64;
%Total points
Mt=m+2;
%h
h=1/(m+1);
%Grid
x=(0:h:1);
y=(0:h:1);

%Intialization of the Matrices
N=Mt*Mt;              %N Unknowns
A=sparse(zeros(N,N)); %Matrix A written in block-diagonal form
F=sparse(zeros(N,1)); %Vector F is the result. The boundary conditions are absorbed into F.

v_exact=sparse(zeros(N,1));
for i=1:Mt
    for j=1:Mt
        n=i+(j-1)*Mt;
        v_exact(n,1)=((exp(pi.*x(i))-exp(-pi.*x(i))).*sin(pi.*y(j)))./(exp(pi)-exp(-pi));
    end
end

for i=1:Mt
    for j=1:Mt
        n=i+(j-1)*Mt;
        G(i,j)=v_exact(n);
    end
end

%Interior Points
for i=2:Mt-1 %Running Loop over x, the first and last grid points (Boundary condtions) are omitted
    for j=2:Mt-1 %Running Loop over y, the first and last grid points (Boundary condtions) are omitted
        n=i+(j-1)*Mt; %Enumeration of the unknown terms in Vector V
        A(n,n)=-4; %Diagonal 
        A(n,n-1)=1; %sub-diagonal
        A(n,n+1)=1; %super-diagonal
        A(n,n-Mt)=1;%far off sub-diagonal
        A(n,n+Mt)=1;%far off super-diagonal
        F(n,1)=0;
    end
end

% Boundary Points: %BC: u(0,y)=0
i=1;
for j=1:Mt
    n=i+(j-1)*Mt;
    A(n,n)=1;
    F(n,1)=0;
end
%BC: u(0,y)=sin
i=Mt;
for j=1:Mt
    n=i+(j-1)*Mt;
    A(n,n)=1;
    F(n,1)=sin(pi*y(j));
end
%BC: u(x,0)=0
j=1;
for i=1:Mt
    n=i+(j-1)*Mt;
    A(n,n)=1;
    F(n,1)=0;
end
%BC: u(x,1)=0
j=Mt;
for i=1:Mt
    n=i+(j-1)*Mt;
    A(n,n)=1;
    F(n,1)=0;
end

V_vec=A\F;

for i=1:Mt
    for j=1:Mt
        n=i+(j-1)*Mt;
        V(i,j)=V_vec(n);
    end
end

x1=(h:h:1-h);
y1=(h:h:1-h);
V1=V(2:Mt-1,2:Mt-1);

figure
mesh(h:h:1-h,h:h:1-h,V1)
title(['Numerical Solution for m = ' num2str(m)'']);
xlabel('x')
ylabel('y')

U1=G(2:Mt-1,2:Mt-1);

figure
mesh(h:h:1-h,h:h:1-h,U1)
title(['Actual solution for m = ' num2str(m)'']);
xlabel('x')
ylabel('y')

figure
mesh(h:h:1-h,h:h:1-h,U1-V1)
title(['Error for m = ' num2str(m)'']);
xlabel('x')
ylabel('y')

Norm_Value=Norm2(U1-V1,h)



function[Norm2value]=Norm2(B,h)
Norm2value=sum(sum(B.^(2)));
Norm2value=h*sqrt(Norm2value);
end



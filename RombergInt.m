function [R,quad,err,h]=RombergInt(f,a,b,n,tol)
%Input    - f is the integrand input as a string 'f'
%         - a and b are uppper and lower limits of integration
%         - n is the maximum number of rows in the table
%         - tol is the tolerance
%Output   - R is the Romberg table
%         - quad is the quadrature value
%         - err is the error estimate
%         - h is the smallest step size used
M=1;
h=b-a;
err=1;
J=0;
R=zeros(4,4);
R(1,1)=h*(feval(f,a)+feval(f,b))/2;
while((err>tol)&&(J<n))||(J<4)
    J=J+1;
    h=h/2;
    s=0;
    for p=1:M
        x=a+h*(2*p-1);
        s=s+feval(f,x);
    end
    R(J+1,1)=R(J,1)/2+h*s;
    M=2*M;
    for K=1:J
        R(J+1,K+1)=R(J+1,K)+(R(J+1,K)-R(J,K))/(4^K-1);
    end
    err=abs(R(J,J)-R(J+1,J+1));
end
quad=R(J+1,J+1);
end
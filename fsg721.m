function dx=fsg721(x)

% 2nd order polynomial
C=[0.107143,0.071429,0.035714];

B=zeros(1,7);
for i=1:3 
    B(i)=C(i); 
end	
B(4)=0.0;
for i=5:7
    B(i)=-C(8-i);
end
A=[1,0];

s=size(x,2);
dx=filter(B,A,x);
dx=[dx(7),dx(7),dx(7),dx(7:s),dx(s),dx(s),dx(s)];

end
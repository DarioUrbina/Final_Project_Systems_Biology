
function dydt = Function_A(t,y,params)

A = y;
%A = inv(A)
%B = y(2);


k0 = params(1,1:11); 
k1f = params(2,1:11);
k1r = params(3,1:11); 
k2 = params(4,1:11); 

%dydt(1,1) = k0 - k2*A - k1f*A + k1r*B;
dydt(1:11,1) = k0 + k2*A - k1f*A + k1r*A/2;

%dydt(2,1) = k1f*A - k1r*B;


end
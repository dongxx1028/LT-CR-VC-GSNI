clear 
clc
% Parameter value
M=8;
n=10^10;
vf1l=5.5;
vf1r=[2.4,4.2,6,7.8,9.6];
v1b=5.8;
vbe=5.6;
zetaL=2;
zetar=2.5;
H=3;
a1=120;
b1=1500;
c1=0.05;
h1=150;
e1=2000;
a2=100;
b2=2000;
c2=-0.08;
h2=200;
e2=1500;
D =5.2*10^(-2);
Lf1r=0.2;
Lf1l=0.25;
L=3;
Ae=600;
Be=500;
Ab=10;
Bb=3;


for i=1:250
    t(i)=0.001*10^(9*i/99);
end
LL = length(t);

% Code for  v(j)(equation (38))
V = zeros(1,M);
for j=1:M
   
    V(j)=0;
    for k=fix((j+1)/2):1:min(j,M/2)
                  V(j)=V(j) + k^(M/2) * factorial (2*k+1) / (factorial(M/2-k+1) * factorial(k+1) * factorial(k) * factorial(j-k+1) * factorial(2*k-j+1));

    end
    V(j)=(-1)^(M/2+j)*V(j);
end


% When vfr=2.4,4.2,6,7.8,9.6, reverse solution of Laplace space  concentration of contaminants at the wellhead (equation (91)) to real space concentration of contaminants
% at the wellhead (equation (92))
for k=1:5
   for i=1:LL
      ft(i,k)=0;
      u=log(2)/t(i);
      for j=1:M
     s=j*u;
% Code for the characteristic roots of the characteristic equation of the left wing of the fracture.: equations£¨34£©-£¨35£©
       lambda1L=(vf1l+sqrt(vf1l^2+4*D*(s-zetaL)))/(2*D);
       lambda2L=(vf1l-sqrt(vf1l^2+4*D*(s-zetaL)))/(2*D);
% Code for the characteristic roots of the characteristic equation of the right wing of the fracture: equations£¨48£©-£¨49£©
       lambda1r=(vf1r(k)+sqrt(vf1r(k)^2+4*D*(s-zetar)))/(2*D);
       lambda2r=(vf1r(k)-sqrt(vf1r(k)^2+4*D*(s-zetar)))/(2*D);
% Code for the characteristic roots of the characteristic equation from the first crack to the bottom of the vertical wellbore: equations£¨68£©-£¨69£© 
       lambda11b=(v1b+sqrt(v1b^2+4*D*s))/(2*D);
       lambda21b=(v1b-sqrt(v1b^2+4*D*s))/(2*D);
% Code for the characteristic roots of the characteristic equation from the bottom of the vertical wellbore to the wellheade: equations£¨73£©-£¨74£© 
       lambda1be=(vbe+sqrt(vbe^2+4*D*s))/(2*D);
       lambda2be=(vbe-sqrt(vbe^2+4*D*s))/(2*D);
% Code for equations£¨80£©-£¨81£©
       d1=a1*(1/(2*n)-1)*Lf1l+e1-c1*h1^2;
       d2=a2*(1-1/(2*n))*Lf1r+e2-c2*h2^2;
% Code for equations£¨82£©-£¨83£©
        q1=c1*(2/(s^3)-2*h1/(s^2)+h1^2/s)+d1/s;
        q2=c2*(2/(s^3)-2*h2/(s^2)+h2^2/s)+d2/s;
% Code for the concentration of contaminanton on the left wing of the first fracture: equation£¨86£©
        Cf1L=(exp(lambda2L*(1-1/(2*n))*Lf1l)*(q1+(lambda1L*a1*(1/(2*n)-1)*Lf1l+a1+e1*lambda1L)/...
           ((lambda2L-lambda1L)*lambda1L^2)+(lambda2L*a1*(1/(2*n)-1)*Lf1l+a1+e1*lambda2L)/((lambda1L-lambda2L)*lambda2L^2))...
           -(a1+e1*lambda1L)/((lambda2L-lambda1L)*lambda1L^2)-(a1+e1*lambda2L)/((lambda1L-lambda2L)*lambda2L^2));
% Code for the concentration of contaminanton the right wing of the first fracture: equation£¨87£©
        Cf1r=(exp(lambda2r*(1/(2*n)-1)*Lf1r)*(q2+(lambda1r*a2*(1-1/(2*n))*Lf1r+a2+e2*lambda1r)/...
           ((lambda2r-lambda1r)*lambda1r^2)+(lambda2r*a2*(1-1/(2*n))*Lf1r+a2+e2*lambda2r)/((lambda1r-lambda2r)*lambda2r^2))...
           -(a2+e2*lambda1r)/((lambda2r-lambda1r)*lambda1r^2)-(a2+e2*lambda2r)/((lambda1r-lambda2r)*lambda2r^2));
% Code for the concentration of contaminanton at the bottom of the wellboree: equation£¨90£©
       C1b=exp(lambda21b*L)*((Cf1L+ Cf1r)/2+(Ab+Bb*lambda11b)/((lambda21b-lambda11b)*lambda11b^2)+(Ab+Bb*lambda21b)/((lambda11b-lambda21b)*lambda21b^2))...
           -(Ab*lambda11b*L+Ab+Bb*lambda11b)/((lambda21b-lambda11b)*lambda11b^2)-(Ab*lambda21b*L+Ab+Bb*lambda21b)/((lambda11b-lambda21b)*lambda21b^2);
% Code for the concentration of contaminanton at the wellhead: equation£¨91£©
       Cbe=exp(lambda2be*H)*(C1b+(Ae+Be*lambda1be)/((lambda2be-lambda1be)*lambda1be^2)+(Ae+Be*lambda2be)/((lambda1be-lambda2be)*lambda2be^2))...
           -(Ae*lambda1be*H+Ae+Be*lambda1be)/((lambda2be-lambda1be)*lambda1be^2)-(Ae*lambda2be*L+Ae+Be*lambda2be)/((lambda1be-lambda2be)*lambda2be^2);
% Code for equation (92)        
       ft(i,k) =  ft(i,k) + V(j) *(Cbe);
      end
     ft(i,k)=ft(i,k)*u;
   end
end
% Double logarithmic curves ofconcentration of contaminanton at the wellhead
loglog(t,ft(:,1),t,ft(:,2),t,ft(:,3),t,ft(:,4),t,ft(:,5), 'k-','Linewidth',1.25)
hold on 
axis([10^(-1) 10^(3) 10^3 10^(5)]);
xlabel('$$t (h)$$ ','FontName','Times New Roman','FontSize',12,'Interpreter','latex') 
ylabel('$$C_{be} (mg/L)$$' ,'FontName','Times New Roman','FontSize',12,'Interpreter','latex')
legend('$$ \overline{v}_{f1r}=2.4 $$','$$ \overline{v}_{f1r}=4.2 $$','$$ \overline{v}_{f1r}=6 $$','$$ \overline{v}_{f1r}=7.8 $$','$$ \overline{v}_{f1r}=9.6 $$','Interpreter','latex')
box off  
grid on   
 
 
 
 
 
 
 
 
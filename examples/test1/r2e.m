%%%%%%%%%%%%%%%%%%%% output: a,e,Omega,i,omega,f,M %%%%%%%%%%%%%%%
function [coe,flag]=r2e(rv,mu)

% rv: position + velocity vector
% mu: standard gravitational parameter

dt=1e-7;
r=rv(1:3)';
v=rv(4:6)';
h=cross(r,v);
i=acos(h(3)/norm(h));

ev=cross(v,h)/mu-r/norm(r);

e=norm(ev);

a=norm(h)^2/(mu*(1-e^2));
n=(sqrt(mu/abs(a)^3));%Æ½½ÇËÙ¶È
if(abs(sin(i))<1e-12)
    omg=0;
else
    omg=atan(-h(1)/h(2));
end

if((h(2)>0))
    omg=omg+pi;
end
   
if(omg<0)
    omg=omg+2*pi;
end



if((abs(sin(i))<1e-12)&&(e>1e-12))
    w=atan(ev(2)/ev(1));u=atan(r(2)/r(1));
    if(ev(1)<0)
        w=w+pi;
    end
    if(r(1)<0)
        u=u+pi;
    end
    
    
elseif((abs(sin(i))>1e-12)&&(e<1e-12))
         w=0;
         u=atan(r(3)/((r(2)*sin(omg)+r(1)*cos(omg))*sin(i)));
        if(r(3)*u<0)
        u=u+pi;
        elseif((r(3)<0)&&(u<0))
        u=u+2*pi;
        end
    
    
elseif((abs(sin(i))<1e-12)&&(e<1e-12))
    w=0;u=atan(r(2)/r(1));
        if(r(1)<0)
        u=u+pi;
        end
    
elseif((abs(sin(i))>1e-12)&&(e>1e-12))
    w=atan(ev(3)/(ev(2)*sin(omg)+ev(1)*cos(omg))/sin(i));
    u=atan(r(3)/((r(2)*sin(omg)+r(1)*cos(omg))*sin(i)));
        if(ev(3)*w<0)
        w=w+pi;
        elseif((ev(3)<0)&&(w<0))
        w=w+2*pi;
        end
    
        if(r(3)*u<0)
        u=u+pi;
        elseif((r(3)<0)&&(u<0))
        u=u+2*pi;
        end
        
        
end

if(e<1)
flag=0;
f=u-w;
f=mod(f,2*pi);


% E=2*atan(sqrt((1-e)/(1+e))*tan(f/2));
% M=E2M(E,e);
% M1=M+n*time;
% E1=M2E(M1,e);
% f1=E2f(E1,e);
% f2=rem(f1+20*pi,2*pi);
% tp=-(E-e*sin(E))/n;
% tp=mod(tp,2*pi/n);

coe=[a e omg i w f];
elseif(e>1)
    flag=1;
    f=acos(   (-a*(e^2-1)/norm(r)-1)/e  );
    rp=r+v*dt;
    if(norm(rp)<norm(r))
        f=-f;
    end
        H=2*atanh(sqrt((e-1)/(1+e))*tan(f/2));
        tp=-(e*sinh(H)-H)/n;
     coe=[a e omg i w f tp];
     
end



function de=f2E(xsit,ecc)
if xsit>=0
    k=fix(0.25*xsit/pi);
else
    k=floor(0.25*xsit/pi);  
end
m=0.5*xsit-2*k*pi; %  m=0.5*xsit+2*k*pi;% end  

if m>=0 && m<0.5*pi    
    de=2*acos(sqrt((((ecc+cos(xsit))/(1+ecc*cos(xsit)))+1)*0.5))+4*k*pi;
elseif  m>=0.5*pi && m<pi    
    de=2*acos(-sqrt((((ecc+cos(xsit))/(1+ecc*cos(xsit)))+1)*0.5))+4*k*pi;
elseif m>=pi && m<1.5*pi    
    de=4*pi-2*acos(-sqrt((((ecc+cos(xsit))/(1+ecc*cos(xsit)))+1)*0.5))+4*k*pi;
elseif m>=1.5*pi && m<2*pi    
    de=4*pi-2*acos(sqrt((((ecc+cos(xsit))/(1+ecc*cos(xsit)))+1)*0.5))+4*k*pi;
end


function M=E2M(E,e)
M=E-e*sin(E);


function E=M2E(M,e)
E=M;
fe=1;
while abs(fe)>1e-15
fe=E-e*sin(E)-M;
dfe=1-e*cos(E);
E=E-fe/dfe;
end


function true=E2f(dde,ecc)
if dde>=0
    k=fix(0.25*dde/pi);
else
    k=floor(0.25*dde/pi);
end
m=0.5*dde-2*k*pi; 

if m>=0 && m<0.5*pi    
    true=2*acos(sqrt((((cos(dde)-ecc)/(1-ecc*cos(dde)))+1)*0.5))+4*k*pi;
elseif  m>=0.5*pi && m<pi    
    true=2*acos(-sqrt((((cos(dde)-ecc)/(1-ecc*cos(dde)))+1)*0.5))+4*k*pi;
elseif m>=pi && m<1.5*pi    
    true=4*pi-2*acos(-sqrt((((cos(dde)-ecc)/(1-ecc*cos(dde)))+1)*0.5))+4*k*pi;
elseif m>=1.5*pi && m<2*pi    
    true=4*pi-2*acos(sqrt((((cos(dde)-ecc)/(1-ecc*cos(dde)))+1)*0.5))+4*k*pi;
end




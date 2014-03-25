% 2d Reaction Diffusion Solver
p.tf=100; p.ht=0.1;p.tt=0:p.ht:p.tf; %
p.xf=2*pi;p.hx=0.4;p.xx=0:p.hx:p.xf;
options=odeset('reltol',1e-8,'abstol',1e-8,'maxstep',1e-1);
p.nt=length(p.tt);p.nx=length(p.xx);
p.alpha=0;p.mu=0;p.beta=1.0;p.w=1.0*ones(p.nx);
p.D=-2*eye(p.nx)+diag(ones(p.nx-1,1),1)+diag(ones(p.nx-1,1),-1);
%Boundary Conditions
%p.D(1,1)=0.0;p.D(1,2)=0.0;p.D(p.nx,p.nx)=0.0;p.D(p.nx,p.nx-1)=0.0;%dirich
%p.D(1,2)=2.0;p.D(p.nx,p.nx-1)=2.0; %reflecting
p.D(1,p.nx)=1;p.D(p.nx,1)=1; %periodic
p.D=.001*p.D/p.hx^2; %diffusivity
u0=(sin((p.xx'*p.xx)/(2*pi))); v0=zeros(p.nx); %data on square
% Set graphics parameters.
fig = figure; set(fig,'color','w') %
p.side=p.nx; x = ((-p.side+1):(p.side))'/p.side;
L=ones(1*length(x)); h = surf(x,x,L); %
[a,e] = view; view(a,90);view(90,90); %
axis([-1 1 -1 1 -2 2]); caxis(26.9*[-1.5 1]);
colormap(hot); shading flat; axis off%
n=0; u=u0; v=v0;

while ishandle(fig)
	u=u+p.ht*(v);
	v=v+p.ht*(p.w-v-(p.beta*cos(p.mu*n*p.ht)).*sin(u)-p.alpha*v.*(v.^2-1)+p.D*u+u*p.D);
	n=n+1;
	if(mod(n,2)==0)
		uu=[cos(u) cos(u); cos(u) cos(u)]; vv=uu; %2x2 patches clarify beh at bndry
		set(h,'zdata',uu,'cdata',10*vv);
	drawnow;
	end;
	if(n>=p.nt)break;end;
end;


% visualization of displacements and von Mises stresses in 2D
% can handle both P1 and P2 elements
% INPUTS:
% p - [2 x Np] matrix of mesh node coordinates
% t - [nbf x Ne] mesh connectivity matrix
% u_h - [Np x 1] vector of mesh node displacements
function plotElasticity(p,t,u_h,D)
t=t'; tFE=t; p=p';
nbf = size(tFE,2);
if nbf==3 % linear finite elements
    pResolution=1;
    shapef=@(r,s)P1shapes(r,s);
elseif nbf==6; % quadratic finite elements
    pResolution=5;
    shapef=@(r,s)P2shapes(r,s);
end
[pbary,tbary]=triBaryGrid(pResolution);

offset=0; T=zeros(0,3); P=zeros(0,2); C=zeros(0,1);
for k=1:size(t,1)
    
    % extracting coordinates of element nodes
    pElmDoFs=t(k,:); xnod=p(pElmDoFs,1); ynod=p(pElmDoFs,2);
    Xnod=p(pElmDoFs,:);
    % extracting displacements at element nodes
    loc2glb=zeros(2*nbf,1);  % local -> global map for degrees of freedom
    loc2glb(1:2:end-1) = 2*t(k,:)-1; loc2glb(2:2:end) = 2*t(k,:);
    unod=u_h(loc2glb);
    
    nsubvert=size(pbary,1);
    x=zeros(nsubvert,1); y=x; z=x; dx=x; dy=x; vmSt=x;
    for i=1:nsubvert
        l1=pbary(i,1); l2=pbary(i,2); l3=pbary(i,3); r=l2; s=l3; 
        [S,dSdr,dSds]=shapef(r,s);
        x(i)=S'*xnod; y(i)=S'*ynod; % global coordinate for (r,s)
        J = [Xnod'*dSdr, Xnod'*dSds]; % Jacobian
        dSdX = J'\[dSdr';dSds'];
        dSdx = dSdX(1,:)'; dSdy = dSdX(2,:)';
        [Phi,B] = basisDisplacementAndStrain(S',dSdx',dSdy');
        d = Phi*unod; % displacement
        dx(i)=d(1); dy(i)=d(2);
        St = D*B*unod; % stress
        vmSt(i)=sqrt(St(1).^2-St(1).*St(2)+St(2).^2+3*St(3).^2);
    end
    T=[T;tbary+offset]; P=[P;[x+dx,y+dy]]; C=[C;vmSt];
    offset=offset+nsubvert;
end
patch('Faces',t(:,1:3),'Vertices',p,'FaceColor','none','EdgeColor',[0.7 0.7 0.7]) % plot initial mesh
patch('Faces',T,'Vertices',P,'FaceVertexCData',C,'FaceColor','interp','EdgeColor','none');
colormap('jet'); colorbar;
% extract and plot triangle edges
hold on
TR=triangulation(T,P); E=freeBoundary(TR)';
x=P(:,1);y=P(:,2);
plot(x(E),y(E),'k-','linewidth',0.5);    
hold off; axis equal;
end


function [Phi,B]=basisDisplacementAndStrain(S,dSdx,dSdy)
nbf=length(S);
Phi=zeros(2,2*nbf); Phi(1,1:2:end)=S; Phi(2,2:2:end)=S;
B=zeros(3,2*nbf); B(1,1:2:end)=dSdx; B(2,2:2:end)=dSdy;
B(3,1:2:end)=dSdy; B(3,2:2:end)=dSdx;
end


% creates a triangle grid of a triangle with grid coordinates expressed
% in barycentric coordinates, n is the number of divisions of a triangle side
function [pbary,tbary]=triBaryGrid(n)
pbary=zeros(0.5*(n+1)*(n+2),3);
% generate barycentric points
c=0;
for i=0:n
    for j=0:n-i
        c = c + 1;
        pbary(c,:) = [(n-i-j)/n,j/n,i/n];
    end
end
% generate triangles
tbary=[];
c=0;
for i=0:n-1
    c=c+1;
    tbary(end+1,:)=[c,c+1,c+1+(n-i)];
    for j=1:n-i-1
        c=c+1;
        tbary(end+1,:)=[c,c+1+(n-i),c+(n-i)];
        tbary(end+1,:)=[c,c+1,c+1+(n-i)];
    end
    c=c+1;
end
end
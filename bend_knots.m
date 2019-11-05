%%%%%%%
% Code to simulate bend knots in circular rope

% Plotting tool obtained from Matlab file exchange, used under the following license:
% 
% Copyright (c) 2016, Janus H. Wesenberg
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% * Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
% 
% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


clear all;

%load new colormaps
colormap(fiber); %load colormap obtained from optomechanical fibers

input = csvread('input_bend_knots.txt', 0, 0, [0,0,0,5]);
h = input(2); %radius of the rod
Eb = input(3)*10^(9);  %bending modulus, moments of inertia real value if E=900kPa
E = input(4)*10^9; %elastic modulus
nu = input(5);
Fpull = input(6)*10^(6); %pulling force



%Other parameters which can be changed

N=10000; %number of timesteps
framenumber=100; %number of frames to capture
dt=10^(-4); %size of timesteps
c_d = 10^5;  %this is the damping parameter given by eta*rho (this is large even when rho is almost zero)


%Parameters derived from inputs

mutwist=Eb/(2*(1+nu));  %shear modulus
Kbulk = E/(3*(1-2*nu)); %bulk modulus
Ar=pi*h^2; %cross-sectional area
I=0.25*pi*h^4; %moment of inertia
J=0.5*pi*h^4; %moment of twist

alpha = Eb*I; %bending energy scale
beta = mutwist*J; %twisting energy scale
gamma = E*Ar; %stretching energy scale

eta2 = c_d * Ar;  %second derivative damping
eta = eta2; %collision stick factor
eta_tw = c_d * J; %twist stick factor



bendList = [
"reef",      "thief",     "granny",    "grief"...
"bowline13", "bowline14", "bowline24", "bowline23"... 
"carrick13", "carrick14"...
"hunter13",  "hunter14",  "hunter24",  "hunter23"...
"zeppelin13","zeppelin14","zeppelin23"...
"ashley13",  "ashley14"...
"alpine13",  "alpine24",  "alpine23"];


bend  = input(1); %input parameter which specifies which knot in knotList to simulate 

load(strcat(pwd, '/initial_conditions/', bendList(bend),'_initial.mat'), 'A', 'B');

fileLabel = strcat(pwd, '/results/', bendList(bend)); %output files will be saved in 'results' subdirectory

A0 = A(:,:,1); B0 = B(:,:,1); %initial quantities from data file

%initial quantities:
n1 = size(A0, 2)-2; n2 = size(B0,2)-2;  n = n1+n2+2; %number of rod elements
LA0=sqrt(sum((A0(:,2:end)-A0(:,1:end-1)).^2)); %unstretched lengths
LB0=sqrt(sum((B0(:,2:end)-B0(:,1:end-1)).^2));
LsA0=Vdom(A0); LsB0=Vdom(B0); %vector containing the mean length of links either side of a dicretized rod element
vel = zeros(3, n+2);

%quantities to track over time
A=zeros(3,n1+2,framenumber); B=zeros(3,n2+2,framenumber); %positions of the knot
TwDispA = zeros(n1,framenumber); TwDispB = zeros(n2, framenumber); %twist displacements



%%%%initialize material frame
d3A = A0(:,2) - A0(:,1);
d1A = cross(cross(d3A, [0,0,1]'), d3A); d1A=d1A./norm(d1A);
d1A = pTransport(A0, d1A);

d3B = B0(:,2) - B0(:,1);
d1B = cross(cross(d3B, [0,0,1]'), d3B); d1B=d1B./norm(d1B);
d1B = pTransport(B0, d1B);

%initial delta twists, and twist angles
mA0 = zeros(1,n1);
thetaA = cumsum([0, mA0]); %initial angle distribution
m_twA = tAngle(A0, d1A, thetaA);
thdotA = zeros(1,n1+1);

mB0 = zeros(1,n2);
thetaB = cumsum([0, mB0]); %initial angle distribution
m_twB = tAngle(B0, d1B, thetaB);
thdotB = zeros(1,n2+1);

Y = A0; Z = B0;
% main loop
for i=1:N 

    %initialize variables, timestepping taking Y0 -> Y, etc.
    Y0 = Y; Z0 = Z; thetaA0 = thetaA; thetaB0 = thetaB;
    LsA = Vdom(Y0); LsB = Vdom(Z0); %vector containing the mean length of links either side of a dicretized rod element

    %contact handling terms
    R = dists([Y0,Z0]); R(:,n1+2) = 100*h; R(n1+2,:) = 100*h;
    iM = interactMat(R,h); pMat = pressure([Y0,Z0],R,iM,h,Kbulk);

    %calculate elastic forces
    FelastA = bForce1(Y0,alpha) + sForce(Y0,LA0,gamma) + tForce(Y0, m_twA, beta);
    FelastB = bForce1(Z0,alpha) + sForce(Z0,LB0,gamma) + tForce(Z0, m_twB, beta);
    iF = iForceB(Y0,Z0,pMat,R); %contact force
    Ftot = [FelastA, FelastB] + iF; %total elastic force
    Ftot(:,1) = [0,0,0]'; Ftot(:,n1+3) = [Ftot(1,n1+3)+Fpull,0,0]'; %enforce BCs


    %calculate hydrodynamic interaction terms implicitly in terms of
    %velocity:
    %calculate interaction strength based on depth of penetration
    Rv = distsV([Y0, Z0]);
    iMc = padarray(iM, [1,1], 'post') + padarray(iM,[1,1],'pre') - padarray(iM,[1,1],'post').*padarray(iM,[1,1],'pre');
    pDepth = Rv.*iMc; 
    pDepth = (pDepth - diag(sum(pDepth))).*[LsA, LsB]';
    Lapl0 = diag([1/LA0(1), 1./LA0(1:end-1)+1./LA0(2:end), 1./LA0(end), 1./LB0(1), 1./LB0(1:end-1)+1./LB0(2:end), 1./LB0(end)]) - diag([0, 1./LA0(2:end), 0, 1./LB0(1), 1./LB0(2:end)], 1) - diag([1./LA0, 0, 1./LB0], -1);     
    Lapl = (eta2 * Lapl0 - eta *pDepth); 
    vel = lsqminnorm(Lapl, Ftot')'; %solve for velocity

    %timestep the positions of the ropes
    Y = Y0 + dt*vel(:,1:n1+2); 
    Z = Z0 + dt*vel(:,n1+3:end);

    %calculate torques
    tauA = angleForce(Y0, m_twA, beta); tauB = angleForce(Z0, m_twB, beta);
    tauA(1) = 0; tauB(1) = 0; %clamped BCs
    tauInt = eta_tw*tauForce([Y0,Z0], vel, iM, R).*[LA0,0,LB0]; %interaction torque term

    %solve implicit equation for angular velocities for rope A
    tauA_tot = tauA + [0, tauInt(2:n1+1)];
    twLA = eta_tw * TwLapl(Y0, LsA0);
    thdotA = lsqminnorm(twLA, tauA_tot')';   
    thetaA = thetaA0 + dt*thdotA;
    d1A = fupdate(Y0, Y, d1A); %update reference frame along the rope
    m_twA = tAngle(Y, d1A, thetaA); %new twist displacement for A

    %solve implicit equation for angular velocities for rope B        
    tauB_tot = tauB + [0, tauInt(n1+4:end)];
    twLB = eta_tw * TwLapl(Z0, LsB0);
    thdotB = lsqminnorm(twLB, tauB_tot')'; 
    thetaB = thetaB0 + dt*thdotB;
    d1B = fupdate(Z0, Z, d1B); %update reference frame along the rope
    m_twB = tAngle(Z, d1B, thetaB); %new twist displacement for B

    %save data every few timesteps
    if mod(i,N/framenumber)==0
        A(:,1:size(Y,2),framenumber*i/N) = Y;    B(:,1:size(Z,2),framenumber*i/N) = Z;
        TwDispA(:,framenumber*i/N) = m_twA';    TwDispB(:, framenumber*i/N) = m_twB';
        fprintf('\nRunning %d%%\n', 100*i/N);
    end
end
% end main loop

%save data
save(strcat(fileLabel,'.mat'), 'A', 'B', 'TwDispA', 'TwDispB');


%%%%%create AVI object
nFrames = size(A,3);
vidObj = VideoWriter(char(strcat(fileLabel,'.avi')));
vidObj.Quality = 100;
vidObj.FrameRate = 20;
open(vidObj);

%%%%create movie
for j=1:1:nFrames

    %Create strain vectors for ropes A and B
    SA = strain(A(:,:,j), TwDispA(:,j)', LA0, h);
    SB = strain(B(:,:,j), TwDispB(:,j)', LB0, h);
    SA = log(1+SA);    SB = log(1+SB);
    SA(1)=SA(2); SA(end)=SA(end-1);
    SB(1)=SB(2); SB(end)=SB(end-1);

    %plot rope A colored by strain
    [x,y,z,C]=tubeplot2(A(:,:,j),h,20,h/5,[0,SA,0]);
    Aplot=surf(x,y,z,C);
    shading interp;
    set(Aplot,'meshstyle','row');
    set(Aplot, 'FaceLighting','phong','SpecularColorReflectance', 0, 'SpecularExponent', 50, 'DiffuseStrength', 1);
    axis([-6,6,-2.5,4,-5,5]);
    view(2); caxis([0,0.7]); c1 = colorbar; ylabel(c1, "Strain");

    %set(gcf, 'Position', get(0,'Screensize')); % Maximize figure
    set(gca, 'XTick',[]); set(gca, 'YTick',[]);

    daspect([1,1,1]); material shiny; camlight right; camlight left;
    grid off; box on;
    set(gca,'Color',[0.6 0.6 0.6]);
    drawnow;
    hold on;

    %plot rope B colored by strain
    [x,y,z,C]=tubeplot2(B(:,:,j),h,20,h/2,[0,SB,0]);
    Bplot=surf(x,y,z,C);
    shading interp;
    set(Bplot,'meshstyle','row');
    set(Bplot, 'FaceLighting','phong','SpecularColorReflectance', 0, 'SpecularExponent', 50, 'DiffuseStrength', 1);
    material shiny;

    drawnow;




    writeVideo(vidObj, getframe(gcf));
    hold off;
    fprintf('\nCompiling video %d%%\n', 100*j/nFrames);
end
close(vidObj);
hold off;






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function R = dists(X) %distance between i'th and j'th links
    Xm = 0.5*(X(:,1:end-1)+X(:,2:end)); %n+1 points, midpoint of each link
    L2 = diag(Xm'*Xm);
    R = -2*(Xm'*Xm) + (L2 + L2');
    R = R.^0.5; 
end

function R = distsV(X) %inverse squared distance between i'th and j'th points
    L2 = diag(X'*X);
    R = -2*(X'*X) + (L2 + L2');
    R = 1./(R + eye(size(R))) - eye(size(R));
end

function iM = interactMat(R, h) %Set R(i,j)=0 for |i-j|<=10
    n = size(R)-1;
    R([1:n+2:end, n+1+1:n+2:end-1, 2*(n+1)+1:n+2:end-2, 3*(n+1)+1:n+2:end-3, 4*(n+1)+1:n+2:end-4, 5*(n+1)+1:n+2:end-5, 6*(n+1)+1:n+2:end-6, 7*(n+1)+1:n+2:end-7, 8*(n+1)+1:n+2:end-8, 9*(n+1)+1:n+2:end-9, 10*(n+1)+1:n+2:end-10])=2*h; 
    R([2:n+2:end-(n+1), 3:n+2:end-2*(n+1), 4:n+2:end-3*(n+1), 5:n+2:end-4*(n+1), 6:n+2:end-5*(n+1), 7:n+2:end-6*(n+1), 8:n+2:end-7*(n+1), 9:n+2:end-8*(n+1), 10:n+2:end-9*(n+1), 11:n+2:end-10*(n+1)])=2*h; 
    iM = R<2*h; %iM(i,j) = 1 iff link i and j are interacting
end


function pMat = pressure(X,R,iM,h,Kbulk)
    R(iM==0) = 0;
    pMat = Kbulk*( 1 - R./(2*h) + 100*(1-R./(2*h)).^3 ); %set cubic term to around 100
    d3 = X(:,2:end) - X(:,1:end-1);
    L = abs(d3'*d3); L(iM==0)=0;
    pMat = pMat.*L;
end


function Fi = iForceB(Y, Z, pMat, R)
    X = [Y,Z]; n1 = size(Y,2)-2;
    Xm = 0.5*(X(:,1:end-1)+X(:,2:end)); %n+1 points, midpoint of each link
    R(R==0) = 1;
    pMat = pMat./R;
    pMat = sum(pMat).*Xm - (Xm*pMat);
    pMat=[pMat(:,1), 0.5*(pMat(:,2:n1+1)+pMat(:,1:n1)), pMat(:,n1+2), pMat(:,n1+3), 0.5*(pMat(:,n1+4:end)+pMat(:,n1+3:end-1)), pMat(:,end)];
    Fi=pMat;
end


function L = TwLapl(X, Ls)
    ph = phi(X); ph = cos(ph(2:end-1));
    L = diag([1./Ls(2), 1./Ls(2:end-2)+1./Ls(3:end-1),1./Ls(end-1)]) - diag([0,ph(2:end)./Ls(3:end-1)],1) - diag(ph(1:end)./Ls(2:end-1),-1);
end


function tau = tauForce(X, F, iM, R)
    d3 = X(:,2:end) - X(:,1:end-1); d3 = d3./vecnorm(d3);
    iM = iM./(eye(size(R)) + R.^3);
    iM = iM .* vecnorm(d3)'; d3 = d3./vecnorm(d3);
    Xm = 0.5*(X(:,1:end-1)+X(:,2:end)); %n+1 points, midpoint of each link
    Fm = 0.5*(F(:,1:end-1)+F(:,2:end));
    tau = cross(Xm, Fm)*iM - cross(Xm, Fm*iM);
    tau = tau - cross(Xm*iM, Fm) + full(sum(iM)).*(cross(Xm, Fm));
    tau = dot(tau, d3);
end


function B = pTransport(X,v0) %X is the set of positions, 3 by n+2, %v0 is initial vector
    n=size(X);
    n=n(2)-2;
    Y=zeros(3,n+1);
    Y(:,1)=v0; %initial frame
    for i=2:n+1
        ti=X(:,i)-X(:,i-1);
        tf=X(:,i+1)-X(:,i);
        R = rot(ti,tf);
        Y(:,i)=R*Y(:,i-1);
    end
    B=Y;
end



function R = rot(ti,tf)
        theta=atan2(norm(cross(ti,tf)),dot(ti,tf)); %angle to rotate ti into tf
        if norm(cross(ti,tf))==0
            m=[0,0,0];
        else
            m=cross(ti,tf)./norm(cross(ti,tf)); %the sign of m comes from here
        end
        m=theta*m;
        A=[0,-m(3),m(2);m(3),0,-m(1);-m(2),m(1),0]; %infinitesimal rotation
        R=expm(A);
end



function kb = darboux(X)
    d3=(X(:,2:end)-X(:,1:end-1));
    t=d3./sqrt(dot(d3,d3));
    kb=2*cross(t(:,1:end-1),t(:,2:end))./(1+dot(t(:,1:end-1),t(:,2:end)));
end



function Fs = sForce(X0,L0,gamma)
    d3 = X0(:,2:end)-X0(:,1:end-1);
    Ln = vecnorm(d3); d3 = d3./Ln;
    stretch = Ln./L0 - 1;
    F = gamma.*(stretch + 0*stretch.^3).*d3;
    Fs = [F, [0,0,0]'] + [[0,0,0]', -F];
end





function ph = phi(X) %turning angle from one segment to the next 
    n=size(X,2)-2;
    d3=X(:,2:n+2)-X(:,1:n+1);
    Y=atan2(sqrt(sum(cross(d3(:,1:n),d3(:,2:n+1)).^2,1)),dot(d3(:,1:n),d3(:,2:n+1)));
    ph=[0,Y,0];
end



function Ls = Vdom(X) %vector containing size of voronoi domain at each X
    n=size(X,2)-2;
    Y=(X-[X(:,1),X(:,1:n+1)]);
    Z=([X(:,2:n+2),X(:,n+2)]-X);
    Ls=0.5*(sqrt(sum(Y.^2))+sqrt(sum(Z.^2)));
end

function k = kappa(X) %vectorize this
    ph=phi(X);
    k=ph./Vdom(X);
end


        
function dp = dphi(X) %dp(i,j,:) = dphi(i)/dx(j)
    n=size(X);
    n=n(2)-2;
    Yp=zeros(n+2,n+2,3);
    d3=X(:,2:n+2)-X(:,1:n+1);
    cosph = (dot(d3(:,1:n), d3(:,2:n+1))./(sqrt(sum(d3(:,1:n).^2)).*sqrt(sum(d3(:,2:n+1).^2))));
    dcos = -1./sqrt(1-cosph.^2);
    ld = -d3(:,2:n+1)./(sqrt(sum(d3(:,1:n).^2)).*sqrt(sum(d3(:,2:n+1).^2))) + (d3(:,1:n)./(sum(d3(:,1:n).^2))).*cosph;
    ud = d3(:,1:n)./(sqrt(sum(d3(:,1:n).^2)).*sqrt(sum(d3(:,2:n+1).^2))) - (d3(:,2:n+1)./(sum(d3(:,2:n+1).^2))).*cosph;
    d=[[0,0,0]',-dcos.*(ld+ud),[0,0,0]'];
    ld=[dcos.*ld,[0,0,0]']; ud=[[0,0,0]',dcos.*ud];
    Yp(:,:,1)=(diag(d(1,:))+diag(ud(1,:),1)+diag(ld(1,:),-1));
    Yp(:,:,2)=(diag(d(2,:))+diag(ud(2,:),1)+diag(ld(2,:),-1));
    Yp(:,:,3)=(diag(d(3,:))+diag(ud(3,:),1)+diag(ld(3,:),-1));
    Yp=real(Yp);
    Yp(isnan(Yp))=0; Yp(isinf(Yp))=0;
    dp=Yp;
end


function derivs = dL(X) %derivs(i,j,:) = dL(i)/dx(j) where L(i) is twice the voronoi domain of i'th vertex
    n=size(X);
    n=n(2)-2;
    Yp=zeros(n+2,n+2,3);
    ld=(X(:,1:n+1)-X(:,2:n+2))./(sqrt(sum((X(:,1:n+1)-X(:,2:n+2)).^2)));
    ud=-ld;
    d=[[0,0,0]',ud]-[ud,[0,0,0]'];
    Yp(:,:,1)=diag(d(1,:))+diag(ud(1,:),1)+diag(ld(1,:),-1);
    Yp(:,:,2)=diag(d(2,:))+diag(ud(2,:),1)+diag(ld(2,:),-1); 
    Yp(:,:,3)=diag(d(3,:))+diag(ud(3,:),1)+diag(ld(3,:),-1);
    derivs=Yp;
end



function Fb=bForce1(X,alpha)
    n=size(X);
    n=n(2)-2;
    ph=phi(X); dph=dphi(X);
    L=2*Vdom(X); dl=dL(X);
    Yb=zeros(3,n+2);
    Yb(1,:)= -((2*alpha*ph)./L)*dph(:,:,1) + ((alpha*ph.^2)./L.^2)*dl(:,:,1);
    Yb(2,:)= -((2*alpha*ph)./L)*dph(:,:,2) + ((alpha*ph.^2)./L.^2)*dl(:,:,2);
    Yb(3,:)= -((2*alpha*ph)./L)*dph(:,:,3) + ((alpha*ph.^2)./L.^2)*dl(:,:,3);
    Fb=Yb;
end


    
function m = tAngle(X, d1, theta) %d1(i) is reference frame vector, theta(i) is angle on link i 
    d3=X(:,2:end)-X(:,1:end-1);
    N=cross(d3(:,1:end-1),d3(:,2:end));
    N=N./sqrt(dot(N,N));
    N(isnan(N))=0;
    d3=d3./sqrt(dot(d3,d3));
    ph=phi(X);
    ph=ph(2:end-1);
    pTrans=cos(ph).*d1(:,1:end-1)+((1-cos(ph)).*dot(N,d1(:,1:end-1))).*N + (sin(ph).*cross(N,d1(:,1:end-1)));  
    pTrans(isnan(pTrans))=d1(isnan(pTrans));
    frameTwist=atan2(dot(d1(:,2:end),cross(d3(:,2:end),pTrans)),dot(d1(:,2:end),pTrans));
    m=theta(2:end)-theta(1:end-1)+frameTwist;
end



function Ft = tForce(X,m,beta) %force on x_i due to twisting, m is given by tAngle
    Ls=Vdom(X);
    kb=darboux(X);
    Es=sqrt(dot(X(:,2:end)-X(:,1:end-1),X(:,2:end)-X(:,1:end-1))); %Es(i) is distance between x_i and x_i+1
    %forces on x_i due to changing the frame due to changing m at i-1, i+1,
    %i respectively
    Fd=0.5*[[0,0,0]', [0,0,0]', (m(1:end)./Ls(2:end-1)).*(kb(:,1:end)./Es(2:end))]; 
    Fu=0.5*[(m(1:end)./Ls(2:end-1)).*(-kb(:,1:end)./Es(1:end-1)), [0,0,0]',[0,0,0]'];
    Fe=0.5*[[0,0,0]', (m(1:end)./Ls(2:end-1)).*(kb(:,1:end)./Es(1:end-1) - kb(:,1:end)./Es(2:end))  ,[0,0,0]'];
    twist1=-beta*(Fd+Fu+Fe);
    %force on x_i due to changing the lengths of links in the twist energy
    %expression
    L=2*Ls; dl=dL(X); 
    twist2=[((beta*[0,m,0].^2)./L.^2)*dl(:,:,1);((beta*[0,m,0].^2)./L.^2)*dl(:,:,2);((beta*[0,m,0].^2)./L.^2)*dl(:,:,3)];
    Yt=twist1+twist2;
    Ft=Yt;
end


function Fa=angleForce(X,m,beta)
    Ls=Vdom(X);
    Fa= -beta*([0,m./Ls(2:end-1)]-[m./Ls(2:end-1),0]);
end

function d1New = fupdate(X0,X1,d1)
    d30=X0(:,2:end)-X0(:,1:end-1); d31=X1(:,2:end)-X1(:,1:end-1);
    ph=atan2(sqrt(sum(cross(d30,d31).^2,1)),dot(d30,d31));
    N=cross(d30,d31);
    N=N./sqrt(dot(N,N));
    Y=cos(ph).*d1+((1-cos(ph)).*dot(N,d1)).*N + (sin(ph).*cross(N,d1));  
    Y(isnan(Y))=d1(isnan(Y));
    d1New=Y;
end



function s = strain(X, m, L0, h) %sqrt strain squared at the boundary of a cross section
    Yb=kappa(X); %curvature
    Yt=[0,m,0]./Vdom(X); %twisting stress
    L=sqrt(sum((X(:,2:end)-X(:,1:end-1)).^2));
    L  = 0.5*[L(1),L(1:end-1)+L(2:end),L(end)]; %voronoi domain-ify the stretch stress
    Lv = 0.5*[L0(1),L0(1:end-1)+L0(2:end),L0(end)];
    Ys = (L./Lv - 1);
    s = sqrt(0.5*h^2*Yt.^2 + (h*Yb + Ys).^2);
end




function [x,y,z,C]=tubeplot2(curve,r,n,ct,S)
% Usage: same as above but this gives colours the rod according to a 
% scalar field S along the curve.

  if nargin<3 || isempty(n), n=8;
     if nargin<2, error('Give at least curve and radius');
     end
  end
  if size(curve,1)~=3
    error('Malformed curve: should be [3,N]');
  end
  if nargin<4 || isempty(ct)
    ct=0.5*r;
  end

  
  %Collapse points within 0.5 r of each other
  npoints=1;
  for k=2:(size(curve,2)-1)
    if norm(curve(:,k)-curve(:,npoints))>ct
      npoints=npoints+1;
      curve(:,npoints)=curve(:,k);
    end
  end
  %Always include endpoint
  if norm(curve(:,end)-curve(:,npoints))>0
    npoints=npoints+1;
    curve(:,npoints)=curve(:,end);
  end

  %deltavecs: average for internal points.
  %           first strecth for endpoitns.
  dv=curve(:,[2:end,end])-curve(:,[1,1:end-1]);

  %make nvec not parallel to dv(:,1)
  nvec=zeros(3,1);
  [~,idx]=min(abs(dv(:,1))); nvec(idx)=1;

  xyz=zeros(3,n+1,npoints+2);
  Col=zeros(3,n+1,npoints+2); 

  %precalculate cos and sing factors:
  cfact=repmat(cos(linspace(0,2*pi,n+1)),[3,1]);
  sfact=repmat(sin(linspace(0,2*pi,n+1)),[3,1]);
  
  %Main loop: propagate the normal (nvec) along the tube
  for k=1:npoints
    convec=cross(nvec,dv(:,k));
    convec=convec./norm(convec);
    nvec=cross(dv(:,k),convec);
    nvec=nvec./norm(nvec);
    %update xyz:
    xyz(:,:,k+1)=repmat(curve(:,k),[1,n+1])+...
        cfact.*repmat(r*nvec,[1,n+1])...
        +sfact.*repmat(r*convec,[1,n+1]);
    Col(:,:,k+1)=S(k);
  end
  %finally, cap the ends:
  xyz(:,:,1)=repmat(curve(:,1),[1,n+1]);
  xyz(:,:,end)=repmat(curve(:,end),[1,n+1]);
  %,extract results:
  x=squeeze(xyz(1,:,:));
  y=squeeze(xyz(2,:,:));
  z=squeeze(xyz(3,:,:));
  Ct=squeeze(Col(1,:,:));
  C=Ct;
  %... and plot:
  if nargout<3, surf(x,y,z,C); end
end






%define realspace position vectors
basvec=[kron(ones(1,Ly),[0:sqrt(3)/2:sqrt(3)*(Lx-1)/2])',(kron(ones(1,Ly/2),[kron(ones(1,Lx/2),[1/2,0]),kron(ones(1,Lx/2),[0,1/2])])+kron([0:3/2:3*(Ly-1)/2],ones(1,Lx)))'];

%find nn and nnn indices without double counting.
%find nn indices:
count=0;
nnindices=[];
for p=1:length(basvec)
    for q=p:length(basvec)
        if abs(sum(abs(basvec(p,:)-basvec(q,:)).^2,2)-1)<0.001
            count=count+1;
            nnindices(count,:)=[p,q];
        end   
    end
end
%PBC nn indices
count=0;
PBCynnindices=[];
for p=1:length(basvec)
    for q=p:length(basvec)
        if (abs(abs(basvec(p,2)-basvec(q,2))-(3*(Ly/2-1)+2))<0.001) && abs(basvec(p,1)-basvec(q,1))<0.001
            count=count+1;
            PBCynnindices(count,:)=[p,q];
        end   
    end
end
count=0;
PBCxnnindices=[];
for p=1:length(basvec)
    for q=p:length(basvec)
        if (abs(abs(basvec(p,1)-basvec(q,1))-(sqrt(3)*(Lx/2-1)+sqrt(3)/2))<0.001) && abs(abs(basvec(p,2)-basvec(q,2))-1/2)<0.001
            count=count+1;
            PBCxnnindices(count,:)=[p,q];
        end   
    end
end
%find nnn indices
count=0;
nnnindices=[];
for p=1:length(basvec)
    for q=p:length(basvec)
        if abs(sum(abs(basvec(p,:)-basvec(q,:)).^2,2)-3)<0.001
            count=count+1;
            nnnindices(count,:)=[p,q];
        end   
    end
end
%find PBC nnn indices
count=0;
PBCxnnnindices=[];
for p=1:length(basvec)
    for q=p:length(basvec)
        if (abs(abs(basvec(p,1)-basvec(q,1))-(sqrt(3)*(Lx/2-1)))<0.001) && abs(abs(basvec(p,2)-basvec(q,2)))<0.001
            count=count+1;
            PBCxnnnindices(count,:)=[p,q];
        elseif (abs(abs(basvec(p,1)-basvec(q,1))-(sqrt(3)*(Lx/2-1)+sqrt(3)/2))<0.001) && abs(abs(basvec(p,2)-basvec(q,2))-3/2)<0.001
            count=count+1;
            PBCxnnnindices(count,:)=[p,q];
        end   
    end
end
count=0;
PBCynnnindices=[];
for p=1:length(basvec)
    for q=p:length(basvec)
        if (abs(abs(basvec(p,1)-basvec(q,1))-(sqrt(3)/2))<0.001) && abs(abs(basvec(p,2)-basvec(q,2))-(3*(Ly/2-1)+3/2))<0.001
            count=count+1;
            PBCynnnindices(count,:)=[p,q];
        end   
    end
end
count=0;
PBCxynnnindices=[];
for p=1:length(basvec)
    for q=p:length(basvec)
        if (abs(abs(basvec(p,1)-basvec(q,1))-(sqrt(3)*(Lx/2-1)+sqrt(3)/2))<0.001) && abs(abs(basvec(p,2)-basvec(q,2))-(3*(Ly/2-1)+3/2))<0.001
            count=count+1;
            PBCxynnnindices(count,:)=[p,q];
        end   
    end
end
%define matrices with the sign
sgnPBCxnnindices=[];
sgnPBCynnindices=[];
sgnPBCxnnnindices=[];
sgnPBCynnnindices=[];
sgnPBCxynnnindices=[];
count=0;
for pp=1:length(PBCxnnindices)
    count=count+1;
    sgnPBCxnnindices(count,1)=sign(basvec(PBCxnnindices(pp,2),1)-basvec(PBCxnnindices(pp,1),1));
    sgnPBCxnnindices(count,2)=sign(basvec(PBCxnnindices(pp,2),2)-basvec(PBCxnnindices(pp,1),2));
end
count=0;
for pp=1:length(PBCynnindices)
    count=count+1;
    sgnPBCynnindices(count,1)=sign(basvec(PBCynnindices(pp,2),1)-basvec(PBCynnindices(pp,1),1));
    sgnPBCynnindices(count,2)=sign(basvec(PBCynnindices(pp,2),2)-basvec(PBCynnindices(pp,1),2));
end
count=0;
for pp=1:length(PBCxnnnindices)
    count=count+1;
    sgnPBCxnnnindices(count,1)=sign(basvec(PBCxnnnindices(pp,2),1)-basvec(PBCxnnnindices(pp,1),1));
    sgnPBCxnnnindices(count,2)=sign(basvec(PBCxnnnindices(pp,2),2)-basvec(PBCxnnnindices(pp,1),2));
end
count=0;
for pp=1:length(PBCynnnindices)
    count=count+1;
    sgnPBCynnnindices(count,1)=sign(basvec(PBCynnnindices(pp,2),1)-basvec(PBCynnnindices(pp,1),1));
    sgnPBCynnnindices(count,2)=sign(basvec(PBCynnnindices(pp,2),2)-basvec(PBCynnnindices(pp,1),2));
end
count=0;
for pp=1:length(PBCxynnnindices)
    count=count+1;
    sgnPBCxynnnindices(count,1)=sign(basvec(PBCxynnnindices(pp,2),1)-basvec(PBCxynnnindices(pp,1),1));
    sgnPBCxynnnindices(count,2)=sign(basvec(PBCxynnnindices(pp,2),2)-basvec(PBCxynnnindices(pp,1),2));
end

%%%%%%%%%%%%%%%%%start loops
for Mchoice=1:length(Marr)
M=Marr(Mchoice);

for Vrandchoice=1:length(Vrandarr)
Vrand=Vrandarr(Vrandchoice)
tnn_dis=tnn_disarr(Vrandchoice);
tnnn_dis=tnnn_disarr(Vrandchoice);
%%% define filenames:
datestring=datestr(now,'yymmddHHMMSS')


%Definitions necessary for bott index 
N1=Lx/2;
N2=Ly/2;
X1=(repmat(rem(1:Lx,2).*[2:Lx+1]/2+((-rem(1:Lx,2)+1).*([1:Lx]/2+0.5)),1,Ly))-1;
Y1temp=[repmat([0.5:3:3*N2;0:3:3*N2-1],N1,1);repmat([1.5:3:3*N2-1;2:3:3*N2],N1,1)];%repmat([2*ceil((1:Ly)/2)+ceil((2:Ly+1)/2-1);1*ceil((1:Ly)/2)+2*ceil((2:Ly+1)/2-1)+0.5],N1,1);
%Y1temp=repmat([0:Ly-1],Lx,1);
Y1=Y1temp(:);
expBIxmat=diag(exp(j*(2*pi/(N1))*X1));
%expBIymat=diag(exp(j*(2*pi/(2*N2))*Y1));
expBIymat=diag(exp(j*(2*pi/(3*N2))*Y1));


Hmass=diag(repmat([M*(kron(ones(1,Lx/2),[1,-1])),M*(kron(ones(1,Lx/2),[-1,1]))],1,Ly/2));
for disavg=1:disavmax
    %disavg
    expH=eye(Lx*Ly);
    Hrand=diag(Vrand*(-0.5+rand(Lx*Ly,1)));
    %define array of nn and nnn couplings
    tnnarr=tnn+tnn_dis*(-0.5+rand(length(nnindices),1));
    tnnPBCxarr=tnn+tnn_dis*(-0.5+rand(length(PBCxnnindices),1));
    tnnPBCyarr=tnn+tnn_dis*(-0.5+rand(length(PBCynnindices),1));
    tnnnarr=tnnn+tnnn_dis*(-0.5+rand(length(nnnindices),1));
    tnnnPBCxarr=tnnn+tnnn_dis*(-0.5+rand(length(PBCxnnnindices),1));
    tnnnPBCyarr=tnnn+tnnn_dis*(-0.5+rand(length(PBCynnnindices),1));
    tnnnPBCxyarr=tnnn+tnnn_dis*(-0.5+rand(length(PBCxynnnindices),1));
    %define Hamiltonian without any magnetic field
    
    
    
    
    for tchoice=1:Tdiv
        tchoice;
        Ax=A*(sin(w*((tchoice)*dt)));
        Ay=A*(cos(w*((tchoice)*dt)));
        H1=zeros(Lx*Ly);
        %nn elements
        for pp=1:length(nnindices)
            H1(nnindices(pp,1),nnindices(pp,2))=tnnarr(pp)*exp(1j*(Ax*(basvec(nnindices(pp,2),1)-basvec(nnindices(pp,1),1))+Ay*(basvec(nnindices(pp,2),2)-basvec(nnindices(pp,1),2))));
        end
        for pp=1:length(PBCxnnindices)
            H1(PBCxnnindices(pp,1),PBCxnnindices(pp,2))=PBCx*tnnPBCxarr(pp)*exp(1j*(Ax*((basvec(PBCxnnindices(pp,2),1)-basvec(PBCxnnindices(pp,1),1))-sgnPBCxnnindices(pp,1)*(Lx/2)*sqrt(3))+Ay*(basvec(PBCxnnindices(pp,2),2)-basvec(PBCxnnindices(pp,1),2))));
        end
        for pp=1:length(PBCynnindices)
            H1(PBCynnindices(pp,1),PBCynnindices(pp,2))=PBCy*tnnPBCyarr(pp)*exp(1j*(Ax*((basvec(PBCynnindices(pp,2),1)-basvec(PBCynnindices(pp,1),1)))+Ay*(basvec(PBCynnindices(pp,2),2)-basvec(PBCynnindices(pp,1),2)-sgnPBCynnindices(pp,2)*3*Ly/2)));
        end
        %nnn indices
        for pp=1:length(nnnindices)
            H1(nnnindices(pp,1),nnnindices(pp,2))=tnnnarr(pp)*exp(1j*(Ax*(basvec(nnnindices(pp,2),1)-basvec(nnnindices(pp,1),1))+Ay*(basvec(nnnindices(pp,2),2)-basvec(nnnindices(pp,1),2))));
        end
        for pp=1:length(PBCxnnnindices)
            H1(PBCxnnnindices(pp,1),PBCxnnnindices(pp,2))=PBCx*tnnnPBCxarr(pp)*exp(1j*(Ax*(basvec(PBCxnnnindices(pp,2),1)-basvec(PBCxnnnindices(pp,1),1)-sgnPBCxnnnindices(pp,1)*sqrt(3)*Lx/2)+Ay*(basvec(PBCxnnnindices(pp,2),2)-basvec(PBCxnnnindices(pp,1),2))));
        end
        for pp=1:length(PBCynnnindices)
            H1(PBCynnnindices(pp,1),PBCynnnindices(pp,2))=PBCy*tnnnPBCyarr(pp)*exp(1j*(Ax*(basvec(PBCynnnindices(pp,2),1)-basvec(PBCynnnindices(pp,1),1))+Ay*(basvec(PBCynnnindices(pp,2),2)-basvec(PBCynnnindices(pp,1),2)-sgnPBCynnnindices(pp,2)*3*Ly/2)));
        end
        for pp=1:length(PBCxynnnindices)
            H1(PBCxynnnindices(pp,1),PBCxynnnindices(pp,2))=PBCx*PBCy*tnnnPBCxyarr(pp)*exp(1j*(Ax*(basvec(PBCxynnnindices(pp,2),1)-basvec(PBCxynnnindices(pp,1),1)-sgnPBCxynnnindices(pp,1)*sqrt(3)*Lx/2)+Ay*(basvec(PBCxynnnindices(pp,2),2)-basvec(PBCxynnnindices(pp,1),2)-sgnPBCxynnnindices(pp,2)*3*Ly/2)));
        end
        H=H1+H1'+Hrand+Hmass;
      

        expH=expH*expm(-1j*H*dt);
    end
    %%%% eigenvalues and eigenvectors for the unitary
    [Utemp,dUtemp]=spdiags(expH);
    szU=size(expH);
    U00=spdiags(Utemp,dUtemp,szU(1),szU(2));
    logU=1i/T*logm(full(U00));
    Hamiltonian=0.5*(logU+logU');
    [Hamiltoniantemp,dHamiltoniantemp]=spdiags(Hamiltonian);
    H00=spdiags(Hamiltoniantemp,dHamiltoniantemp,szU(1),szU(2));
    [W,d]=eig(full(H00));
    %d=eig(full(H00));
    En(:,disavg)=diag(d);   
    %bott index         
   for movingboundchoice=1:length(movingboundarr)
       movingbound=movingboundarr(movingboundchoice);
       p2=max(find(diag(d)<=max(movingbound,fixedbound)));
       p1=min(find(diag(d)>min(movingbound,fixedbound)));
      UX=W(:,p1:p2)'*(expBIxmat)*W(:,p1:p2);
       UY=W(:,p1:p2)'*(expBIymat)*W(:,p1:p2);
       Ubott=UY*UX*UY'*UX';
       [Ubotttemp,dUbotttemp]=spdiags(Ubott);
       szUbott=size(Ubott);
       Ubott00=spdiags(Ubotttemp,dUbotttemp,szUbott(1),szUbott(2));
       index(movingboundchoice,disavg)=imag(sum(log(eig(full(Ubott00)))))/(2*pi);
   end
Name=sprintf('data/graphenefloquetdisorderdata%s-%s.mat',datestring,JobID);
save(Name,'Lx','Ly','PBCx','PBCy','A','M','w','tnn','tnnn','tnn_dis','tnnn_dis' ,'T','Tdiv','dt','Vrand','disavmax','seedvalue','En','index','movingboundarr')%,'avgIPR')
   
end
Name=sprintf('data/graphenefloquetdisorderdata%s-%s.mat',datestring,JobID);
save(Name,'Lx','Ly','PBCx','PBCy','A','M','w','tnn','tnnn','tnn_dis','tnnn_dis' ,'T','Tdiv','dt','Vrand','disavmax','seedvalue','En','index','movingboundarr')%,'avgIPR')
toc

end
end

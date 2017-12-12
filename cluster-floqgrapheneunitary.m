for Mchoice=1:length(Marr)
M=Marr(Mchoice)
for Vrandchoice=1:length(Vrandarr)
Vrand=Vrandarr(Vrandchoice)
%%% define filenames:
datestring=datestr(now,'yymmddHHMMSS');
sprintf('%s-%s',datestring,JobID)


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
avgIPR=zeros(Lx*Ly,length(movingboundarr));
for disavg=1:disavmax
    %disavg
    expH=eye(Lx*Ly);
    Hrand=diag(Vrand*(-0.5+rand(Lx*Ly,1)));
    for tchoice=1:Tdiv
        %tchoice
        Ax=A*(sin(w*((tchoice)*dt)));
        Ay=A*(cos(w*((tchoice)*dt)));
        tnn1=tnn*exp(1j*0.5*(Ax*sqrt(3)+Ay));
        tnn2=tnn*exp(1j*0.5*(Ax*sqrt(3)-Ay));
        tnn3=tnn*exp(-1j*Ay);
        tnn1conj=tnn*exp(-1j*0.5*(Ax*sqrt(3)+Ay));
        tnn2conj=tnn*exp(-1j*0.5*(Ax*sqrt(3)-Ay));
        tnn3conj=tnn*exp(1j*Ay);
        tnnn1=tnnn*exp(1j*Ax*sqrt(3));
        tnnn2=tnnn*exp(-1j*0.5*(Ax*sqrt(3)+3*Ay));
        tnnn3=tnnn*exp(1j*0.5*(sqrt(3)*Ax-3*Ay));
        tnnn1conj=tnnn1';
        tnnn2conj=tnnn2';
        tnnn3conj=tnnn3';

        %define relevant block matrices (see note)

        H11=diag([kron(ones(1,Lx/2-1),[tnn1,tnn2]),tnn1],1)+diag(tnnn1*ones(1,Lx-2),2)+PBCx*diag(tnnn1conj*[1,1],Lx-2)+PBCx*diag(tnn2conj,Lx-1);
        H11tot=H11+H11';
        H22=diag([kron(ones(1,Lx/2-1),[tnn2,tnn1]),tnn2],1)+diag(tnnn1*ones(1,Lx-2),2)+PBCx*diag(tnnn1conj*[1,1],Lx-2)+PBCx*diag(tnn1conj,Lx-1);
        H22tot=H22+H22';
        H12=diag(kron(ones(1,Lx/2),[tnn3,0]))+diag(tnnn3*ones(1,Lx-1),1)+diag(tnnn2*ones(1,Lx-1),-1)+PBCx*diag(tnnn2,Lx-1);%+PBCx*diag(tnnn3,-Lx+1);
        H12conj=H12';
        H21=diag(kron(ones(1,Lx/2),[0,tnn3]))+diag(tnnn3*ones(1,Lx-1),1)+diag(tnnn2*ones(1,Lx-1),-1)+PBCx*diag(tnnn2,Lx-1);%+PBCx*diag(tnnn3,-Lx+1);
        H21conj=H21';
        %time-dependent Hamiltoniaan
        H=Hmass+Hrand+kron(diag(kron(ones(1,Ly/2),[1,0])),H11tot)+kron(diag(kron(ones(1,Ly/2),[0,1])),H22tot)+kron(diag([kron(ones(1,Ly/2-1),[1,0]),1],1),H12)+kron(diag([kron(ones(1,Ly/2-1),[0,1]),0],1),H21)+kron(diag([kron(ones(1,Ly/2-1),[1,0]),1],-1),H12conj)+kron(diag([kron(ones(1,Ly/2-1),[0,1]),0],-1),H21conj)+PBCy*kron(diag(1,Ly-1),H21conj)+PBCy*kron(diag(1,-Ly+1),H21);
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
    %%% quasienergy dependent quantities    
    for movingboundchoice=1:length(movingboundarr)
        %bott index
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
        %avgIPR calculated at the different quasienergies
        eigenveclist=find(abs(diag(d)-movingbound)<energywidthtolerance);
        avgIPR(:,movingboundchoice)=avgIPR(:,movingboundchoice)+sum(abs(W(:,eigenveclist)).^4,2)/length(eigenveclist);
    end
    
end
avgIPR=avgIPR/disavmax;
Name=sprintf('data/graphenefloquetdisorderdata%s-%s.mat',datestring,JobID);
save(Name,'Lx','Ly','PBCx','PBCy','A','M','w','tnn','tnnn','T','Tdiv','dt','Vrand','disavmax','seedvalue','En','index','avgIPR','movingboundarr','energywidthtolerance')
toc

end
end
%dos
%figure()
%[e1,h1]=hist(En(:),35)
%plot(h1,e1,'-o')

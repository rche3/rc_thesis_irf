function [bigD,gamma0,vbias1,vbias2,vbias] = proc_vbias(iflag,nztwid,phitwid,omegaUtwid,nq,nk,nN,PX,Eetaztwid,EXeta,delta,bigD,gamma0,vbias1,vbias2)

%local variables used in this function:
%int: kd nz m i j h r s 
%vector[real]: delta0m(nztwid) vec_red_Q(nk^2)
%symmetric[real] Q(nk,nk)
%rect[real] omegaU phi Pxm(nk,nztwid) temp_mat temp_F F_rbtwid F_r F_rtwid
%           vec_Frtwid F_ratwid F_rbtwid vec_Frbtwid Prj Prm red_Q(nk,nk) L_k(nk^2,kd) d_D F gamma0twid temp_B vec_delta vec_gamma0twid vec_gamma0

kd = nk*(nk+1)/2;
nz = nN*nztwid;

%compute phi, omegaU, gamma0twid and gamma0
phi = kron(eye(nN),phitwid);
omegaU = kron(eye(nN),omegaUtwid);

%Note you should calculate gamma0twid and gamma0 in the following steps
%instead of using %psdinit(),otherwise it reports error when phitwid (and phi) is close to identity matrix
vec_gamma0twid = (eye(nztwid^2)-kron(phitwid,phitwid))\omegaUtwid(:);
gamma0twid = reshape(vec_gamma0twid,nztwid,[]);
vec_gamma0 = (eye(nz^2)-kron(phi,phi))\omegaU(:);
gamma0 = reshape(vec_gamma0,nz,[]);

%compute bigD
bigD = PX*gamma0*PX';

%compute vbias1
if iflag==0
	vec_delta = zeros(nz,1);
	for m=0:nq
        temp = delta(:,:,m+1);
		vec_delta = vec_delta + temp(:);
    end
    vbias1 = zeros(nk,1)-((bigD)\PX)*((eye(nz)-phi)\(omegaU*vec_delta));
else
    temp = Eetaztwid(:,:,nq+1)';
	vbias1 = zeros(nk,1)-((bigD)\(EXeta+PX*((eye(nz)-phi)\temp(:))));
end

%compute symmetric matrix Q and redundant Q
%Note for matlab version: in matlab, Q and red_Q are the same
Q=zeros(nk,nk);
for i=1:nk
    for j=1:i
        Q(i,j) = (j-1)*(2*nk-j+2)/2+i-j+1;
        Q(j,i) = Q(i,j);
    end
end
red_Q = Q;

%compute L_k
L_k = zeros(nk^2,kd);
vec_red_Q = red_Q(:);
for i=1:nk^2
	L_k(i,vec_red_Q(i)) = 1;
end

%compute d_D
d_D = zeros(nk^2,kd)-kron(inv(bigD),inv(bigD))*L_k;

%compute F
for r=1:nk
    %compute F_r i.e., F_1, F_2, ... F_k
	if iflag == 0
        vec_Frtwid = zeros(nztwid^2,1);
		for m=1:nN
			Prm = (PX(r:r,(m-1)*nztwid+1:m*nztwid))';

			for s=0:nq
                temp = (phitwid^(s+1))*gamma0twid*Prm*delta(:,m,s+1)'*omegaUtwid;
				vec_Frtwid = vec_Frtwid + temp(:);
            end
        end
		vec_Frtwid = (eye(nztwid^2)-kron(phitwid,phitwid))\vec_Frtwid;
		F_rtwid = reshape(vec_Frtwid,nztwid,[]);
		F_r = kron(eye(nN),F_rtwid);
    else
		%compute F_ratwid
        F_ratwid = zeros(nztwid,nztwid);
		if nq>0
			for j=1:nN
				Prj = PX(r:r,(j-1)*nztwid+1:j*nztwid)';

				for s=0:nq-1
					F_ratwid = F_ratwid + phitwid^(s+1)*gamma0twid*Prj*Eetaztwid(j,:,s+1);
                end
            end
        end

		%compute F_rbtwid
		vec_Frbtwid = zeros(nztwid^2,1);
		for j=1:nN
			Prj = PX(r,(j-1)*nztwid+1:j*nztwid)';
            
            temp = phitwid^(nq+1)*gamma0twid*Prj*Eetaztwid(j,:,nq+1);
			vec_Frbtwid = vec_Frbtwid + temp(:);
        end
		vec_Frbtwid = (eye(nztwid^2)-kron(phitwid,phitwid))\vec_Frbtwid;
		F_rbtwid = reshape(vec_Frbtwid,nztwid,[]);

		%compute F_r i.e., F_1, F_2, ... F_k
		F_r = kron(eye(nN),F_ratwid)+kron(eye(nN),F_rbtwid);
    end

	%compute temp_F, i.e., the r-th column of F
	temp_F = zeros(kd,1);
	for i=1:nk
		for h=i:nk
			temp_F(Q(i,h),1) = PX(i,:)*(F_r+F_r')*PX(h,:)';
        end
    end

	if r==1
		F = temp_F;
	else
		F = [F,temp_F];
    end
end

%compute vbias2 and vbias
for i=1:nk
    temp_B = d_D((i-1)*nk+1:i*nk,1:kd);

	if i==1
		vbias2 = trace(temp_B*F);
	else
		vbias2 = [vbias2;trace(temp_B*F)];
    end
end
vbias = vbias1 + vbias2;
    
end

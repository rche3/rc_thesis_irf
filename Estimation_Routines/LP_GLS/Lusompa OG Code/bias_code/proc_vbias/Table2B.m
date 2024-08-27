clear;
clc;

%1.

%2. Declare and set the parameters.

global nz nq nk nN iflag varphi theta sigma2_1 sigma2_2 phi varphi_vec q_vec qplus1_vec omegaU gamma0 PX bigD delta Eetaz EXeta vbias1 vbias2 vbias Table2B

nz=2;
nk=1;
nN=1;
iflag=0;

varphi_vec=[0.4;-0.09];
q_vec=[0;1;3;7;11;19;39];
qplus1_vec=q_vec+1;

theta=0.3;

sigma2_1=1;
sigma2_2=1;
omegaU=[sigma2_1,0;0,sigma2_2];

EXeta=zeros(nk,1);

%3. Execute the program

Table2B=zeros(size(varphi_vec,1),1+size(q_vec,1));

for i=1:size(varphi_vec,1)
    Table2B(i,1)=varphi_vec(i);
  for j=1:size(q_vec,1)
        varphi=varphi_vec(i);
        phi=[varphi,theta;0,0];
        nq=q_vec(j);

        PX=[0,1];

        delta=zeros(nz,nN,nq+1);
        Eetaz=zeros(nN,nz,nq+1);

        %compute delta
        for h=1:nq+1
            temp=[1,0]*inv(eye(nz)-phi)*(eye(nz)-phi^(nq+2-h));
            delta(:,:,h)=temp';
        end

        vbias = proc_vb_d(nz,phi,omegaU,nk,PX,delta,nq);

        Table2B(i,j+1)=vbias;
    end
end

Table2B

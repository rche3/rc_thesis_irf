source proc_vbias.src
source proc_vb_d.src

*1. State the purpose of the run.

display "This run is in response to vbias39.pdf, and produces Table2C."
display "It was exectuted at" %DATEANDTIME()
display " "


*2. Declare and set the parameters.

declare int nz nq nk nN iflag h i j

compute nz=2
compute nk=2
compute nN=1
compute iflag=0

declare real varphi theta sigma2_1 sigma2_2
declare rect[real] phi omegaU(nz,nz) phi(nz,nz) PX(nk,nz) vbias_matrix
declare vector[real] varphi_vec vbias
declare vector[int] q_vec qplus1_vec
declare vector[rect] delta

compute varphi_vec=|| 0.4 | -0.09 ||
compute q_vec=|| 0 | 1 | 3 | 7 | 11 | 19 | 39 ||
dim qplus1_vec(%rows(q_vec))
do h=1,%rows(q_vec),1
  compute qplus1_vec(h)=q_vec(h)+1
end do h

compute theta=0.3

compute sigma2_1=1
compute sigma2_2=1
compute omegaU=|| sigma2_1,0 | 0,sigma2_2 ||


*3. Execute the program

dim vbias_matrix(nk*%rows(varphi_vec),1+%rows(q_vec))

do i=1,%rows(varphi_vec),1
  do j=1,%rows(q_vec),1
    compute varphi=varphi_vec(i)
    compute phi=|| varphi,theta | 0,0 ||
    compute nq=q_vec(j)

    compute PX=%zeros(nz,nz)
    do h=1,nq+1,1
      compute PX=PX+phi^h
    end do h
    compute PX=(|| 1.0,0 ||*PX)~~|| 0,1.0 ||

    dim delta(nq+1)

    * compute delta
    do h=1,nq+1,1
      compute delta(h)=tr(|| 1.0,0 ||*inv(%identity(nz)-phi)*(%identity(nz)-phi^(nq+2-h)))
    end do h

    @proc_vb_d nZ phi omegaU nk PX delta nq vbias

    compute vbias_matrix(i,1)=varphi
    compute vbias_matrix(i,1+j)=vbias(1)
    compute vbias_matrix(nk+i,1)=varphi
    compute vbias_matrix(nk+i,1+j)=vbias(2)
  end do j
end do i

report(action=define,hlabels=|| "","varphi","","","q+1","","","",""||)
report(action=modify,row=input,atcol=1,fillby=rows) "" "" qplus1_vec
do h=1,nk*%rows(varphi_vec),1
  report(action=modify,row=new,atcol=1,fillby=rows) "vbias" %xrow(vbias_matrix,h)
end do h
report(action=format,picture="*.#",align=right)
report(action=show)


source proc_vbias.src

procedure proc_vb_d nZtwid phitwid omegaUtwid nk PX delta nq vbias

  type int nztwid nq nk
  type rect[real] omegaUtwid phitwid PX
  type vector[real]   *vbias
  type vector[rect] delta

  local vect[real] EXeta vbias1 vbias2
  local vect[rect] Eetaztwid
  local rect gamma0 bigD
  local int nN iflag

  comp nN=1, iflag=0
  dim vbias1(nk) vbias2(nk) gamma0(nztwid,nztwid) bigD(nk,nk)
  dim EXeta(nk) Eetaztwid(nq+1)

  @proc_vbias iflag nztwid phitwid omegaUtwid nq nk nN PX Eetaztwid EXeta delta bigD gamma0 vbias1 vbias2 vbias

  return
end



*1. State the purpose of the run.

display "This run is in response to vbias39.pdf, and produces Table2A."
display "It was exectuted at" %DATEANDTIME()
display " "


*2. Declare and set the parameters.

declare int nz nq nk nN iflag h i j

compute nz=2
compute nq=1
compute nk=1
compute nN=1
compute iflag=0

declare real varphi theta sigma2_1 sigma2_2
declare rect[real] phi omegaU(nz,nz) gamma0 phi(nz,nz) PX(nk,nz) bigD
declare vector[real] varphi_vec EXeta(nk) vbias1 vbias2 vbias vbias_vec vbias1_vec vbias2_vec
declare vector[int] q_vec qplus1_vec
declare vector[rect] delta Eetaz

compute varphi_vec=|| 0.4 | -0.4 | -0.09 ||
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

dim vbias_vec(%rows(q_vec)) vbias1_vec(%rows(q_vec)) vbias2_vec(%rows(q_vec))
report(action=define,hlabels=|| "","varphi","","","q+1","","","",""||)
report(action=modify,row=input,atcol=1,fillby=rows) "" "" qplus1_vec

do i=1,%rows(varphi_vec),1
  do j=1,%rows(q_vec),1
    compute varphi=varphi_vec(i)
    compute phi=|| varphi,theta | 0,0 ||
    compute nq=q_vec(j)

    compute PX=%zeros(nz,nz)
    do h=1,nq+1,1
      compute PX=PX+phi^h
    end do h
    compute PX=|| 1.0,0 ||*PX

    dim delta(nq+1) Eetaz(nq+1)

    * compute delta
    do h=1,nq+1,1
      compute delta(h)=tr(|| 1.0,0 ||*inv(%identity(nz)-phi)*(%identity(nz)-phi^(nq+2-h)))
    end do h

*    @proc_vbias iflag nz phi omegaU nq nk nN PX Eetaz EXeta delta bigD gamma0 vbias1 vbias2 vbias
     @proc_vb_d nZ phi omegaU nk PX delta nq vbias

    compute vbias_vec(j)=%scalar(vbias)
    compute vbias1_vec(j)=%scalar(vbias1)
    compute vbias2_vec(j)=%scalar(vbias2)
  end do j

  report(action=modify,row=new,atcol=1,fillby=rows) "vbias" varphi vbias_vec
*  report(action=modify,row=new,atcol=1,fillby=rows) "vbias1" "" vbias1_vec
*  report(action=modify,row=new,atcol=1,fillby=rows) "vbias2" "" vbias2_vec
  report(action=format,picture="*.##",align=right)
end do i

report(action=show)

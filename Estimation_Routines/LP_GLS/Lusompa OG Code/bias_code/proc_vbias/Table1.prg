source proc_vbias.src
source proc_vb_ma0.src
source proc_vb_d.src

*the code is set up to be extended to produce bias for q>0, for the MA process eta(t)=v(t)+...+v(t+q), v(t)~iid.
*But Table 1 only has results for q=0.

*1. State the purpose of the run.

display "This run is in response to vbias39.pdf, and produces Table1."
display "It was exectuted at" %DATEANDTIME()
display " "


*2. Declare and set the parameters.

declare int nz nq nk nN iflag g h i j

compute nz=2
compute nk=1
compute nN=1
compute iflag=0

declare real phi1A phi2A phi1B phi2B phi1C phi2C sigma2_u
declare rect[real] phi_matrix omegaU(nz,nz) phi(nz,nz) PX(nk,nz) Table1_vb_d Table1_ma0 Eetaztwid0(1,nz)
declare vector[real] vbias
declare vector[int] q_vec
declare vector[rect] delta Eetaz


compute phi1A = 1.5
compute phi2A = -0.54
compute phi1B = 0.55
compute phi2B = 0.4
compute phi1C = 0.9
compute phi2C = 0
compute phi_matrix=|| phi1A, phi2A | phi1B, phi2B | phi1C, phi2C ||

compute q_vec=|| 0 ||

compute sigma2_u=1
compute omegaU=%zeros(nz,nz)
compute omegaU(1,1)=sigma2_u

compute PX=%zeros(nk,nz)
compute PX(1,1)=1


*3. Execute the program

dim Table1_vb_d(%rows(phi_matrix),%rows(q_vec)) Table1_ma0(%rows(phi_matrix),%rows(q_vec))

do g=1,%rows(phi_matrix),1

    compute phi = %xrow(phi_matrix,g)~~|| 1.0,0.0 ||
    compute nq=q_vec(1)

    dim delta(nq+1) Eetaz(nq+1)
    * compute delta
    do h=1,nq+1,1
      do i=1,nz,1
        {
*         defining Eetaz in this way sets delta=1
          if i==1
            compute Eetaz(h) = || sigma2_u ||
          else
            compute Eetaz(h) = Eetaz(h)~||0.0||
        }
        {
          if i==1
            compute delta(h) = || 1.0 ||
          else
            compute delta(h) = delta(h)~~||0.0||
        }

      end do i
    end do h
    comp EetaZtwid0 = Eetaz(1)

    @proc_vb_d nz phi omegaU nk PX delta nq vbias
    compute Table1_vb_d(g,1)=%scalar(vbias)
    @proc_vb_ma0 nz phi omegaU nk PX EetaZtwid0 vbias
    compute Table1_ma0(g,1)=%scalar(vbias)

end do g

display "Table 1B, col (2), using proc_vb_d" Table1_vb_d
display "Table 1B, col (2), using proc_vb_ma0" Table1_ma0


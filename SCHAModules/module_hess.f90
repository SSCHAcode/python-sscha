! This subroutine calculates the stochastic average of an element of the
! third order force constants based on fu2 and lmat. 
! The force constants element is given with three indexes that correspond 
! to atom-cartesian coordinates.
!
! NOTE: there is a first addend that is not implemented at the moment
!       because it should be zero if the odd correction is calculated
!       at the atomic positions that minimize the gradient with respect
!       to Wyckoff positions. This non-implemented addend vanishes
!       if positions are fixed by symmetry.
! Version de get_v3_element_sym pero guardando datos en una matrix auxiliar

module module_hess

use omp_lib
implicit none

contains

subroutine get_v3_red(nat, norbit3, orbit3t, orbit3o, indep_3fc_elem, n_indep_3fc_elem, kernel_3fc, rot_3fc, &
    ur, eprod, f, rho, log_err, s_cart, irt, translations_irt, v3_red, ref_3fc, &
    nref3, dim2, nat_sc, n_mode, n_random, nr, nperm, nsym)

    use stochastic
    implicit none

    integer, intent(in) :: nat
    integer, dimension(nref3), intent(in) :: norbit3, n_indep_3fc_elem
    integer, dimension(nref3,dim2,3), intent(in) :: orbit3t
    integer, dimension(nref3,dim2,2), intent(in) :: orbit3o
    double precision, dimension(nref3,27,27), intent(in) :: kernel_3fc
    double precision, dimension(nperm, nsym, 27, 27), intent(in) :: rot_3fc
    integer, dimension(nref3,27), intent(in) :: indep_3fc_elem

    double precision, dimension(n_random,n_mode), intent(in) :: ur
    double precision, dimension(n_mode,n_mode), intent(in) :: eprod
    double precision, dimension(n_random,nat_sc,3), intent(in) :: f
    double precision, dimension(n_random), intent(in) :: rho
    character (len=10), intent(in) :: log_err
    double precision, dimension(3,3,48), intent(in) :: s_cart
    integer, dimension(48,nat_sc), intent(in) :: irt
    integer, dimension(nat_sc,nr), intent(in) :: translations_irt

    double precision, dimension(nat*3,nat_sc*3,nat_sc*3), intent(out) :: v3_red
    double precision, dimension(nref3, 27), intent(out) :: ref_3fc

    integer :: nref3, dim2, nat_sc, n_mode, n_random, nr, nperm, nsym
    double precision :: v3
    double precision, dimension(maxval(n_indep_3fc_elem)) :: indep_3fc
    double precision, dimension(27) :: tmp_ref_3fc
    double precision, dimension(27) :: aux_3fc

    integer :: ref3, equiv, nat1, nat2, nat3, i, alpha, beta, gamma, ind1, ind2, ind3, comb_ind1, comb_ind2, index, iperm, isym


    indep_3fc = 0 ! We need to initialise it to zero.
    ref_3fc = 0
    tmp_ref_3fc = 0
    v3_red = 0
    aux_3fc = 0

    !$omp parallel private (equiv, i, nat1, nat2, nat3, alpha, beta, gamma, v3, indep_3fc, tmp_ref_3fc)
    !$omp do schedule (dynamic, 1) private (ind1, comb_ind1, ind2, comb_ind2, ind3, index, iperm, isym, aux_3fc)
    do ref3 = 1, nref3
            do equiv = 1, norbit3(ref3)
                    nat1 = orbit3t(ref3,equiv,1)
                    nat2 = orbit3t(ref3,equiv,2)
                    nat3 = orbit3t(ref3,equiv,3)
                    if (equiv == 1) then
                            do i = 1, n_indep_3fc_elem(ref3)
                                    alpha = indep_3fc_elem(ref3,i)/9
                                    beta = mod(indep_3fc_elem(ref3,i),9)/3
                                    gamma = mod(indep_3fc_elem(ref3,i),3)
                                    call get_v3_element_sym(3*nat1+alpha, 3*nat2+beta, 3*nat3+gamma, &
                            ur, eprod, f, rho, log_err, nsym, s_cart, irt, translations_irt, v3, &
                            nat_sc, n_mode, n_random, nr)
                                    indep_3fc(i) = v3
                            call dgemv('N', 27, n_indep_3fc_elem(ref3), 1.0d0, &
                            kernel_3fc(ref3, :, :n_indep_3fc_elem(ref3)), 27, &
                    indep_3fc(:n_indep_3fc_elem(ref3)), 1, 0.0d0, tmp_ref_3fc, 1)
                            ref_3fc(ref3, :) = tmp_ref_3fc
                            end do
                            do ind1 = 0, 2
                                    comb_ind1 = 3*nat1+ind1+1
                                    do ind2 = 0, 2
                                            comb_ind2 = 3*nat2+ind2+1
                                            do ind3 = 0, 2
                                                    index = (3*ind1+ind2)*3+ind3+1
                                                    v3_red(comb_ind1, comb_ind2, 3*nat3+ind3+1) = &
                                            ref_3fc(ref3, index)
                                            end do
                                    end do
                            end do
                    else
                            iperm = orbit3o(ref3, equiv, 1)+1
                            isym = orbit3o(ref3, equiv, 2)+1
                            call dgemv('N', 27, 27, 1.0d0, rot_3fc(iperm,isym,:,:), 27, ref_3fc(ref3,:), 1, &
                    0.0d0, aux_3fc, 1)
                            do ind1 = 0, 2
                                    comb_ind1 = 3*nat1+ind1+1
                                    do ind2 = 0, 2
                                            comb_ind2 = 3*nat2+ind2+1
                                            do ind3 = 0, 2
                                                    index = (3*ind1+ind2)*3+ind3+1
                                                    v3_red(comb_ind1, comb_ind2, 3*nat3+ind3+1) = aux_3fc(index)
                                            end do
                                    end do
                            end do
                    end if
            end do
    end do
    !$omp end do
    !$omp end parallel
  end subroutine get_v3_red

subroutine get_ref3fc(nat, orbit3t, indep_3fc_elem, n_indep_3fc_elem, kernel_3fc, rot_3fc,&
    ur, eprod, f, rho, log_err, s_cart, irt, translations_irt, verbose, ref_3fc, &
    nref3, dim2, nat_sc, n_mode, n_random, nr, nperm, nsym)

    use stochastic
    implicit none

    integer, intent(in) :: nat
    integer, dimension(nref3), intent(in) :: n_indep_3fc_elem
    integer, dimension(nref3,dim2,3), intent(in) :: orbit3t
    double precision, dimension(nref3,27,27), intent(in) :: kernel_3fc
    double precision, dimension(nperm, nsym, 27, 27), intent(in) :: rot_3fc
    integer, dimension(nref3,27), intent(in) :: indep_3fc_elem

    double precision, dimension(n_random,n_mode), intent(in) :: ur
    double precision, dimension(n_mode,n_mode), intent(in) :: eprod
    double precision, dimension(n_random,nat_sc,3), intent(in) :: f
    double precision, dimension(n_random), intent(in) :: rho
    character (len=10), intent(in) :: log_err
    double precision, dimension(3,3,48), intent(in) :: s_cart
    integer, dimension(48,nat_sc), intent(in) :: irt
    integer, dimension(nat_sc,nr), intent(in) :: translations_irt
    
    logical, intent(in) :: verbose

    double precision, dimension(nref3, 27), intent(out) :: ref_3fc

    integer :: nref3, dim2, nat_sc, n_mode, n_random, nr
    double precision :: v3
    double precision, dimension(maxval(n_indep_3fc_elem)) :: indep_3fc
    double precision, dimension(27) :: tmp_ref_3fc

    integer :: ref3, nat1, nat2, nat3, i, alpha, beta, gamma, nperm, nsym


    indep_3fc = 0 ! We need to initialise it to zero.
    ref_3fc = 0
    tmp_ref_3fc = 0

    if (verbose) then
        print *, "======================= get_ref3fc() ======================= "
        print *, "Computing the independent elements for the third order FCs."
    end if

    !$omp parallel private (i, nat1, nat2, nat3, alpha, beta, gamma, v3, indep_3fc, tmp_ref_3fc)
    !$omp do schedule (dynamic, 1) private (v3, indep_3fc, tmp_ref_3fc)
    do ref3 = 1, nref3
        nat1 = orbit3t(ref3,1,1)
        nat2 = orbit3t(ref3,1,2)
        nat3 = orbit3t(ref3,1,3)
        do i = 1, n_indep_3fc_elem(ref3)
                alpha = indep_3fc_elem(ref3,i)/9
                beta = mod(indep_3fc_elem(ref3,i),9)/3
                gamma = mod(indep_3fc_elem(ref3,i),3)
                call get_v3_element_sym(3*nat1+alpha, 3*nat2+beta, 3*nat3+gamma, &
        ur, eprod, f, rho, log_err, nsym, s_cart, irt, translations_irt, v3, &
        nat_sc, n_mode, n_random, nr)
                indep_3fc(i) = v3
        call dgemv('N', 27, n_indep_3fc_elem(ref3), 1.0d0, &
        kernel_3fc(ref3, :, :n_indep_3fc_elem(ref3)), 27, &
indep_3fc(:n_indep_3fc_elem(ref3)), 1, 0.0d0, tmp_ref_3fc, 1)
        ref_3fc(ref3, :) = tmp_ref_3fc
        end do
    end do
    !$omp end do
    !$omp end parallel
    if (verbose) then
        print *, "=======================     DONE     ======================="
        print *, ""
    end if
end subroutine get_ref3fc

subroutine get_ref4fc(orbit4t, indep_4fc_elem, n_indep_4fc_elem, kernel_4fc, rot_4fc, &
    ur, eprod, f, rho, log_err, s_cart, irt, translations_irt, verbose, ref_4fc, &
    nref4, dim2, nat_sc, n_mode, n_random, nr, nperm, nsym)

    use stochastic
    implicit none

    integer, dimension(nref4), intent(in) :: n_indep_4fc_elem
    integer, dimension(nref4,dim2,4), intent(in) :: orbit4t
    double precision, dimension(nref4,81,81), intent(in) :: kernel_4fc
    double precision, dimension(nperm, nsym, 81, 81), intent(in) :: rot_4fc
    integer, dimension(nref4,81), intent(in) :: indep_4fc_elem

    double precision, dimension(n_random,n_mode), intent(in) :: ur
    double precision, dimension(n_mode,n_mode), intent(in) :: eprod
    double precision, dimension(n_random,nat_sc,3), intent(in) :: f
    double precision, dimension(n_random), intent(in) :: rho
    character (len=10), intent(in) :: log_err
    double precision, dimension(3,3,48), intent(in) :: s_cart
    integer, dimension(48,nat_sc), intent(in) :: irt
    integer, dimension(nat_sc,nr), intent(in) :: translations_irt

    logical, intent(in) :: verbose
    double precision, dimension(nref4, 81), intent(out) :: ref_4fc

    integer :: nref4, dim2, nat_sc, n_mode, n_random, nr
    double precision :: v4
    double precision, dimension(maxval(n_indep_4fc_elem)) :: indep_4fc
    double precision, dimension(81) :: tmp_ref_4fc

    integer :: ref4, equiv, nat1, nat2, nat3, nat4, i, alpha, beta, gamma, delta
    integer :: nperm, nsym


    indep_4fc = 0 ! We need to initialise it to zero.
    ref_4fc = 0
    tmp_ref_4fc = 0

    if (verbose) then
        print *, "======================= get_ref4fc() ======================= "
        print *, "  Computing the independent elements for the 4th order FCs"
    end if

    !$omp parallel private (i, nat1, nat2, nat3, nat4, alpha, beta, gamma, delta)
    !$omp do schedule (dynamic, 1) private (v4, indep_4fc, tmp_ref_4fc)
    do ref4 = 1, nref4
        nat1 = orbit4t(ref4,1,1)
        nat2 = orbit4t(ref4,1,2)
        nat3 = orbit4t(ref4,1,3)
        nat4 = orbit4t(ref4,1,4)
        do i = 1, n_indep_4fc_elem(ref4)
            alpha = indep_4fc_elem(ref4,i)/27
            beta = mod(indep_4fc_elem(ref4,i),27)/9
            gamma = mod(mod(indep_4fc_elem(ref4,i),27),9)/3
            delta = mod(indep_4fc_elem(ref4,i),3)
            call get_v4_element_sym(3*nat1+alpha, 3*nat2+beta, 3*nat3+gamma, 3*nat4+delta, &
ur, eprod, f, rho, log_err, nsym, s_cart, irt, translations_irt, v4, &
nat_sc, n_mode, n_random, nr)
            indep_4fc(i) = v4
        call dgemv('N', 81, n_indep_4fc_elem(ref4), 1.0d0, &
        kernel_4fc(ref4, :, :n_indep_4fc_elem(ref4)), 81, &
indep_4fc(:n_indep_4fc_elem(ref4)), 1, 0.0d0, tmp_ref_4fc, 1)
        ref_4fc(ref4, :) = tmp_ref_4fc
        end do
    end do
    !$omp end do
    !$omp end parallel
    if (verbose) then
        print *, "=======================     DONE     ======================="
        print *, ""
    end if
end subroutine get_ref4fc

subroutine get_ref_vsq(refq2,pol_vecs,rot_3fc,ref_3fc,mapping_triplet,verbose,v1,&
	n_mode_sc,n_mode,iq,nperm,nsym,nref3,nrefq2,dimq2,nat,nat_sc)

	implicit none

	integer, dimension(nrefq2,dimq2,2), intent(in) :: refq2

	complex, dimension(iq, n_mode, n_mode_sc), intent(in) :: pol_vecs

	double precision, dimension(nperm, nsym, 27, 27), intent(in) :: rot_3fc
    double precision, dimension(nref3, 27), intent(in) :: ref_3fc

    integer, dimension(nat,nat_sc,nat_sc,3), intent(in) :: mapping_triplet

    logical, intent(in) :: verbose

	complex, dimension(nrefq2, n_mode, n_mode, n_mode), intent(out) :: v1

	double precision, dimension(27) :: aux_3fc
    double precision :: tstart, tend
	integer :: n_mode_sc, n_mode, ns, qmu, qnu, mu, nu, at3, iq, index, nrefq2, dimq2, i,j, rq2
	integer :: nperm,nsym,nref3,ref3,equiv,iperm,isym,a,b,c,alpha,beta,gamma,nat,nat_sc
	complex, dimension(:), allocatable :: laux1, laux2
    complex, dimension(:,:), allocatable :: lres1

	logical, parameter :: debug = .false.

	ns = n_mode_sc

	! Allocate stuff

	allocate(lres1(n_mode,n_mode_sc))
    allocate(laux1(n_mode_sc))
	allocate(laux2(n_mode_sc))

    if (verbose) then
        tstart = omp_get_wtime()
        print*, "======================= get_ref_vsq() ======================"
        print*, ""
        print*, "        Computing all reference V_{a}^{\alpha}(q1,q2)"
        print*, ""
    end if
    v1 = 0
    !$omp parallel private (qmu,qnu,at3,a,alpha,mu,laux1,lres1,ref3,equiv,iperm,isym,i,j)
    !$omp do schedule (dynamic, 1) private (b,c,aux_3fc,beta,gamma,index,nu,laux2)
    do rq2 = 1, nrefq2
        print "(A, I8, A, I8)", &
"    Calculating V(rq2=", rq2,") out of", nrefq2
        qmu = refq2(rq2,1,1)
        qnu = refq2(rq2,1,2)
        do at3 = 1, n_mode
            a = (at3-1)/3
            alpha = mod(at3-1,3)
            lres1(:,:) = (0.0d0, 0.0d0)
            do b = 0, nat_sc-1
                do c = 0, nat_sc-1
                    ref3 = mapping_triplet(a+1,b+1,c+1,1)+1
                    iperm = mapping_triplet(a+1,b+1,c+1,2)+1
                    isym = mapping_triplet(a+1,b+1,c+1,3)+1
                    aux_3fc = 0.0d0
                    do j = 1, 27
                        do i = 1, 27
                            aux_3fc(i) = aux_3fc(i) + &
rot_3fc(iperm, isym, i, j) * ref_3fc(ref3, j)
                        end do
                    end do
                    do mu = 1, n_mode
                        do beta = 0, 2
                            do gamma = 0, 2
                                index = &
(3*alpha+beta)*3+gamma+1
                                lres1(mu,3*c+gamma+1) = &
lres1(mu,3*c+gamma+1) + aux_3fc(index) * pol_vecs(qmu+1,mu,3*b+beta+1)
                            end do
                        end do
                    end do
                end do
            end do
            do mu = 1, n_mode
                do nu = 1, n_mode  
                    laux2 = CONJG(pol_vecs(qnu+1,nu,:))
                    v1(rq2, at3, mu, nu) = sum(lres1(mu,:)*laux2) 
                end do
            end do
        end do
    end do
    !$omp end do
    !$omp end parallel
	if (verbose) then
        tend = omp_get_wtime()
        print*, ""
        print "(A, ES16.6,A)", "  Elapsed time inside get_ref_vsq():", tend-tstart, " seconds"
		print*, "=======================     DONE     ======================="
        print*, ""
	endif
    deallocate(laux1)
	deallocate(lres1)
	deallocate(laux2)

end subroutine get_ref_vsq

subroutine get_all_vsq(pol_vecs, v3_red, map_uc, vs, iq, n_mode, n_mode_sc, nat_sc)

    implicit none

    complex(8), dimension(iq,n_mode,n_mode_sc), intent(in) :: pol_vecs
    double precision, dimension(n_mode,n_mode_sc,n_mode_sc), intent(in) :: v3_red

    integer, dimension(nat_sc, nat_sc), intent(in) :: map_uc

    complex(8), dimension(n_mode_sc,iq,iq,n_mode,n_mode), intent(out) :: vs


    integer :: iq, n_mode, n_mode_sc, nat_sc, nat

    integer :: q1,q2,mu1,mu2,a,b,c,at,bt,ct
    integer :: alpha,beta,gamma

    vs = 0.0d0
    do concurrent(q1=1:iq,q2=1:iq)
        do concurrent(mu1=1:n_mode,mu2=1:n_mode)
            do concurrent(at=1:nat_sc,bt=1:nat_sc,ct=1:nat_sc)
                a = map_uc(at,at)+1
                b = map_uc(at,bt)+1
                c = map_uc(at,ct)+1
                do concurrent(alpha=1:3,beta=1:3,gamma=1:3)
                    vs(3*(at-1)+alpha,q1,q2,mu1,mu2)=vs(3*(at-1)+alpha,q1,q2,mu1,mu2)+&
v3_red(3*(a-1)+alpha,3*(b-1)+beta,3*(c-1)+gamma)*&
pol_vecs(q1,mu1,3*(bt-1)+beta)*conjg(pol_vecs(q2,mu2,3*(ct-1)+gamma))
                end do
            end do
        end do
    end do

end subroutine get_all_vsq

subroutine get_ref_wsq(refq4,pol_vecs,rot_4fc,ref_4fc,mapping_quadruplet,verbose,&
    v1,n_mode_sc,n_mode,iq,nperm,nsym,nref4,nrefq4,dimq4,nat,nat_sc)

	implicit none

	integer, dimension(nrefq4,dimq4,4), intent(in) :: refq4

	complex(8), dimension(iq, n_mode, n_mode_sc), intent(in) :: pol_vecs

	double precision, dimension(nperm, nsym, 81, 81), intent(in) :: rot_4fc
    double precision, dimension(nref4, 81), intent(in) :: ref_4fc

    integer, dimension(nat,nat_sc,nat_sc,nat_sc,3), intent(in) :: mapping_quadruplet
    
    logical, intent(in) :: verbose

	complex(8), dimension(nrefq4, n_mode, n_mode, n_mode, n_mode), intent(out) :: v1

	double precision, dimension(81) :: aux_4fc
	double precision, dimension(3,3,3,3) :: aux_4fc_3

    complex(8), dimension(9,9) :: aux_4fc_ij
	integer :: n_mode_sc, n_mode, ns, q1, q2, q3, q4, iq, nrefq4, dimq4, ii, jj, rq4
	integer :: nperm,nsym,nref4,ref4,iperm,isym,at,bt,ct,dt,alpha,beta,gamma,delta,nat,nat_sc
    integer :: mu1,mu2,mu3,mu4,a,b,c,d

    double precision, parameter :: pi  = 4.0d0*atan(1.0d0)
    double precision :: tstart, tend
    complex(8), parameter :: j = (0.0d0, 1.0d0)
    complex(8), dimension(n_mode*n_mode,n_mode*n_mode) :: vflat
    complex(8), dimension(nat_sc,n_mode*n_mode,n_mode*n_mode) :: vtmp
    complex(8), dimension(9,n_mode*n_mode) :: tmp
    complex(8), dimension(nat_sc*nat_sc,9,n_mode*n_mode) :: L, R

    if (verbose) then
        tstart = omp_get_wtime()
        print*, "======================= get_ref_wsq() ======================"
        print*, ""
        print "(A, I8, A)", &
"   Computing the", nrefq4, " reference W^{(1)}(q1,q2,q3,q4)"
        print*, ""
        print*, " This might take some time..."
        print*, ""
    end if

	ns = n_mode_sc

	v1 = 0.0d0
    tstart = omp_get_wtime()
    do rq4 = 1, nrefq4
        print "(A, I8, A, I8)", "    Calculating W(rq4=", rq4,") out of", nrefq4
		q1 = refq4(rq4,1,1)+1
		q2 = refq4(rq4,1,2)+1
		q3 = refq4(rq4,1,3)+1
		q4 = refq4(rq4,1,4)+1
        !$omp parallel private (bt,alpha,beta)
        !$omp do schedule (static) private (mu1,mu2)
        do at = 1, nat_sc
            do bt=1, nat_sc
                do alpha = 1, 3
                    do mu1 = 1, n_mode
                        do beta = 1, 3
                            do mu2 = 1, n_mode
                                L(nat_sc*(at-1)+bt,3*(alpha-1)+beta,n_mode*(mu1-1)+mu2) = &
CONJG(pol_vecs(q1,mu1,3*(at-1)+alpha))*pol_vecs(q2,mu2,3*(bt-1)+beta)
                            end do
                        end do
                    end do
                end do
            end do
        end do
        !$omp end do
        !$omp end parallel
        !$omp parallel private (dt,gamma,delta)
        !$omp do schedule (static) private (mu3,mu4)
        do ct=1, nat_sc
            do dt=1, nat_sc
                do gamma = 1, 3
                    do mu3 = 1, n_mode
                        do delta = 1, 3
                            do mu4 = 1, n_mode
                                R(nat_sc*(ct-1)+dt,3*(gamma-1)+delta,n_mode*(mu3-1)+mu4) = &
pol_vecs(q3,mu3,3*(ct-1)+gamma)*conjg(pol_vecs(q4,mu4,3*(dt-1)+delta))!*phase
                            end do
                        end do
                    end do
                end do
            end do
        end do
        !$omp end do
        !$omp end parallel
        vflat = 0.0d0
        tmp = 0.0d0
        do a = 1, nat
            vtmp = 0.0d0     
            !$omp parallel private (c,d,tmp,ref4,iperm,isym,aux_4fc,aux_4fc_ij)
            !$omp do schedule (static) private (alpha,beta,gamma,delta)
            do b=1, nat_sc
                do c=1, nat_sc
                    do d=1, nat_sc
                        ref4 = mapping_quadruplet(a,b,c,d,1)+1
                        iperm = mapping_quadruplet(a,b,c,d,2)+1
                        isym = mapping_quadruplet(a,b,c,d,3)+1
                        aux_4fc=0.0d0
                        call dgemv( &
'N',81,81,1.0d0,rot_4fc(iperm,isym,:,:),81,ref_4fc(ref4,:),1,0.0d0,aux_4fc,1)
                        do alpha= 1, 3
                            do beta = 1, 3
                                do gamma = 1, 3
                                        do delta = 1, 3
                                        aux_4fc_ij(3*(alpha-1)+beta,3*(gamma-1)+delta) = &
aux_4fc(3*(3*(3*(alpha-1)+beta-1)+gamma-1)+delta)
                                        end do
                                end do
                            end do
                        end do
                        call zgemm( &
'N','N', 9, n_mode*n_mode, 9, (1.0d0,0.0d0), aux_4fc_ij, 9, R(nat_sc*(c-1)+d,:,:), 9, (0.0d0,0.0d0), tmp, 9)
                        call zgemm( &
'T','N', n_mode*n_mode, n_mode*n_mode, 9, (1.0d0,0.0d0), L(nat_sc*(a-1)+b,:,:), 9, tmp, 9, &
(1.0d0,0.0d0), vtmp(b,:,:), n_mode*n_mode)
                    end do
                end do
            end do
            !$omp end do
            !$omp end parallel
            vflat = vflat + sum(vtmp, dim=1)*iq
        end do
        do mu1=1,n_mode
            do mu2=1,n_mode
                do mu3=1,n_mode
                    do mu4=1,n_mode
                        v1(rq4,mu1,mu2,mu3,mu4) = vflat(n_mode*(mu1-1)+mu2,n_mode*(mu3-1)+mu4)
                    end do
                end do
            end do
        end do
	end do
    if (verbose) then
        tend = omp_get_wtime()
        print*, ""
        print "(A, ES16.6,A)", "  Elapsed time inside get_ref_wsq():", tend-tstart, " seconds"
		print*, "=======================     DONE     ======================="
        print*, ""
	endif
end subroutine get_ref_wsq

subroutine get_scf_wsq(wsq1, F, refq4, refq4o, norbitq4, P, degs, verbose, eps, alpha_mix, &
wsq_scf, nrefq4, n_mode, iq, dimq4, nsym)
    
    implicit none

    complex(8), dimension(nrefq4,n_mode,n_mode,n_mode,n_mode), intent(in) :: wsq1

    double precision, dimension(iq, iq, n_mode, n_mode), intent(in) :: F

	integer, dimension(nrefq4,dimq4,4), intent(in) :: refq4
    integer, dimension(nrefq4,dimq4,2), intent(in) :: refq4o
    integer, dimension(nrefq4), intent(in) :: norbitq4

	complex(8), dimension(iq, n_mode, n_mode, nsym), intent(in) :: P
    logical, dimension(iq, n_mode, n_mode), intent(in) :: degs
    logical, intent(in) :: verbose

    double precision, intent(in) :: eps, alpha_mix

    complex(8), dimension(nrefq4,n_mode,n_mode,n_mode,n_mode), intent(out) :: wsq_scf

    complex(8), dimension(nrefq4,n_mode,n_mode,n_mode,n_mode) :: wsq2
	complex(8), dimension(nrefq4,n_mode,n_mode,n_mode,n_mode) :: wsqn

    double precision :: tstart, tend

    integer, dimension(8,4) :: perms
    integer, dimension(4) :: nunu, lala, mumu1, mumu2
    integer :: nrefq4, n_mode, iq, dimq4, nsym
    integer :: mu1,mu2,mu3,mu4,mu5,mu6,iter
    integer :: rq4_0,rq4_1,rq4_2,q4i_1,q4i_2,q1,q2,q3,q4,q5_1,q6_1,q5_2,q6_2
	integer :: q1_ref, q2_ref, q3_ref, q4_ref, q5_1ref, q5_2ref, q6_1ref, q6_2ref
	integer :: nu1,nu2,nu5,nu6,lambda3,lambda4,lambda5,lambda6
    integer :: q10,q20,q30,q40,q5,q6,iperm_1,iperm_2,isym_1,isym_2

    logical :: its_conv
	complex(8) :: ktea1, ktea2, W1, W2
    integer, parameter :: maxiter = 10

    perms = reshape([ &
        1,2,3,4, 2,1,4,3, 3,4,1,2, 4,3,2,1, &
        1,3,2,4, 3,1,4,2, 2,4,1,3, 4,2,3,1], &
    [8,4], order=[2,1])

    wsq_scf = wsq1

    its_conv = .true.
    if (verbose) then
        tstart = omp_get_wtime()
        print*, "======================= get_scf_wsq() ======================="
        print*, ""
        print*, "        Self-consintent loop to compute W(q1,q2,q3,q4)"
        print*, ""
        print "(A, ES16.6, A, ES16.6)", &
"        eps=", eps, ";  alpha_mix=", alpha_mix
        print*, ""
    end if
    do1 : do iter = 1, maxiter
        if (verbose) then
            print "(A, I4)", "     Iteration", iter
            print*, "    ---------------------"
        end if
        if (iter == 1) then 
            wsq2 = wsq1
        else
            wsq2 = wsq_scf
        end if
        wsqn = (0.0d0,0.0d0)
        do rq4_0 = 1, nrefq4
            q10 = refq4(rq4_0,1,1)+1
            q20 = refq4(rq4_0,1,2)+1
            q30 = refq4(rq4_0,1,3)+1
            q40 = refq4(rq4_0,1,4)+1
            do rq4_1 = 1, nrefq4
                q1_ref = refq4(rq4_1,1,1)+1
                q2_ref = refq4(rq4_1,1,2)+1
                q5_1ref = refq4(rq4_1,1,3)+1
                q6_1ref = refq4(rq4_1,1,4)+1
                do q4i_1 = 1, norbitq4(rq4_1)
                    q1 = refq4(rq4_1,q4i_1,1)+1
                    q2 = refq4(rq4_1,q4i_1,2)+1
                    q5_1 = refq4(rq4_1,q4i_1,3)+1
                    q6_1 = refq4(rq4_1,q4i_1,4)+1
                    iperm_1 = refq4o(rq4_1,q4i_1,1)
					isym_1 = refq4o(rq4_1,q4i_1,2)+1
                    if ((q10 == q1) .and. (q20 == q2)) then
                        do rq4_2 = 1, nrefq4	
                            q5_2ref = refq4(rq4_2,1,1)+1
                            q6_2ref = refq4(rq4_2,1,2)+1
                            q3_ref = refq4(rq4_2,1,3)+1
                            q4_ref = refq4(rq4_2,1,4)+1	
                            do q4i_2 = 1, norbitq4(rq4_2)
                                q5_2 = refq4(rq4_2,q4i_2,1)+1
                                q6_2 = refq4(rq4_2,q4i_2,2)+1
                                q3 = refq4(rq4_2,q4i_2,3)+1
                                q4 = refq4(rq4_2,q4i_2,4)+1
                                iperm_2 = refq4o(rq4_2,q4i_2,1)
								isym_2 = refq4o(rq4_2,q4i_2,2)+1
                                if ((q30 == q3) .and. (q40 == q4)) then
                                    if ((q5_1 == q5_2) .and. (q6_1 == q6_2)) then
                                        !$omp parallel &
                                        !$omp private (mu2,mu3,mu4,mu5,mu6) &
                                        !$omp private (nu1,nu2,nu5,nu6,lambda5,lambda6,lambda3,lambda4,nunu,lala,mumu1,mumu2)
                                        !$omp do schedule (dynamic, 1) private (ktea1,ktea2,W1,W2)         
                                        do mu1 = 1, n_mode
                                        do nu1 = 1, n_mode
                                        if (degs(q1,nu1,mu1)) then
                                        do mu2 = 1, n_mode
                                        do nu2 = 1, n_mode
                                        if (degs(q2,nu2,mu2)) then
                                        do mu5 = 1, n_mode
                                        do nu5 = 1, n_mode
                                        if (degs(q5_1,nu5,mu5)) then
                                        do mu6 = 1, n_mode
                                        mumu1 = [mu1,mu2,mu5,mu6]
                                        do nu6 = 1, n_mode
                                        nunu = [nu1,nu2,nu5,nu6]
                                        if (degs(q6_1,nu6,mu6)) then
                                            ktea1 = &
P(q1_ref,mumu1(perms(iperm_1+1,1)),nunu(perms(iperm_1+1,1)),isym_1)*&
CONJG(P(q2_ref,mumu1(perms(iperm_1+1,2)),nunu(perms(iperm_1+1,2)),isym_1))*&
CONJG(P(q5_1ref,mumu1(perms(iperm_1+1,3)),nunu(perms(iperm_1+1,3)),isym_1))*&
P(q6_1ref,mumu1(perms(iperm_1+1,4)),nunu(perms(iperm_1+1,4)),isym_1)
                                            if (iperm_1 == 0) then
                                                W1 = &
ktea1*wsq2(rq4_1,nunu(perms(iperm_1+1,1)),nunu(perms(iperm_1+1,2)),nunu(perms(iperm_1+1,3)),nunu(perms(iperm_1+1,4)))
                                            elseif (iperm_1 == 1) then
                                                W1 = &
CONJG(ktea1*wsq2(rq4_1,nunu(perms(iperm_1+1,1)),nunu(perms(iperm_1+1,2)),nunu(perms(iperm_1+1,3)),nunu(perms(iperm_1+1,4))))
                                            elseif (iperm_1 == 2) then
                                                W1 = &
CONJG(ktea1*wsq2(rq4_1,nunu(perms(iperm_1+1,1)),nunu(perms(iperm_1+1,2)),nunu(perms(iperm_1+1,3)),nunu(perms(iperm_1+1,4))))
                                            elseif (iperm_1 == 3) then
                                                W1 = &
ktea1*wsq2(rq4_1,nunu(perms(iperm_1+1,1)),nunu(perms(iperm_1+1,2)),nunu(perms(iperm_1+1,3)),nunu(perms(iperm_1+1,4)))
                                            elseif (iperm_1 == 4) then
                                                W1 = &
ktea1*wsq2(rq4_1,nunu(perms(iperm_1+1,1)),nunu(perms(iperm_1+1,2)),nunu(perms(iperm_1+1,3)),nunu(perms(iperm_1+1,4)))
                                            elseif (iperm_1 == 5) then
                                                W1 = &
CONJG(ktea1*wsq2(rq4_1,nunu(perms(iperm_1+1,1)),nunu(perms(iperm_1+1,2)),nunu(perms(iperm_1+1,3)),nunu(perms(iperm_1+1,4))))
                                            elseif (iperm_1 == 6) then
                                                W1 = &
CONJG(ktea1*wsq2(rq4_1,nunu(perms(iperm_1+1,1)),nunu(perms(iperm_1+1,2)),nunu(perms(iperm_1+1,3)),nunu(perms(iperm_1+1,4))))
                                            elseif (iperm_1 == 7) then
                                                W1 = &
ktea1*wsq2(rq4_1,nunu(perms(iperm_1+1,1)),nunu(perms(iperm_1+1,2)),nunu(perms(iperm_1+1,3)),nunu(perms(iperm_1+1,4)))
                                            end if
                                            do mu3 = 1, n_mode
                                            do lambda3 = 1, n_mode
                                            if (degs(q3,lambda3,mu3)) then
                                            do mu4 = 1, n_mode
                                            do lambda4 = 1, n_mode
                                            mumu2 = [mu5,mu6,mu3,mu4]
                                            if (degs(q4,lambda4,mu4)) then
                                            do lambda5 = 1, n_mode
                                            if (degs(q5_2,lambda5,mu5)) then
                                            do lambda6 = 1, n_mode
                                            lala = [lambda5,lambda6,lambda3,lambda4]
                                            if (degs(q6_2,lambda6,mu6)) then
                                                ktea2 = &
P(q5_2ref,mumu2(perms(iperm_2+1,1)),lala(perms(iperm_2+1,1)),isym_2)*&
CONJG(P(q6_2ref,mumu2(perms(iperm_2+1,2)),lala(perms(iperm_2+1,2)),isym_2))*&
CONJG(P(q3_ref,mumu2(perms(iperm_2+1,3)),lala(perms(iperm_2+1,3)),isym_2))*&
P(q4_ref,mumu2(perms(iperm_2+1,4)),lala(perms(iperm_2+1,4)),isym_2)
                                                if (iperm_2 == 0) then
                                                    W2 = &
ktea2*wsq1(rq4_2,lala(perms(iperm_2+1,1)),lala(perms(iperm_2+1,2)),lala(perms(iperm_2+1,3)),lala(perms(iperm_2+1,4)))
                                                elseif (iperm_2 == 1) then
                                                    W2 = &
CONJG(ktea2*wsq1(rq4_2,lala(perms(iperm_2+1,1)),lala(perms(iperm_2+1,2)),lala(perms(iperm_2+1,3)),lala(perms(iperm_2+1,4))))
                                                elseif (iperm_2 == 2) then
                                                    W2 = &
CONJG(ktea2*wsq1(rq4_2,lala(perms(iperm_2+1,1)),lala(perms(iperm_2+1,2)),lala(perms(iperm_2+1,3)),lala(perms(iperm_2+1,4))))
                                                elseif (iperm_2 == 3) then
                                                    W2 = &
ktea2*wsq1(rq4_2,lala(perms(iperm_2+1,1)),lala(perms(iperm_2+1,2)),lala(perms(iperm_2+1,3)),lala(perms(iperm_2+1,4)))
                                                elseif (iperm_2 == 4) then
                                                    W2 = &
ktea2*wsq1(rq4_2,lala(perms(iperm_2+1,1)),lala(perms(iperm_2+1,2)),lala(perms(iperm_2+1,3)),lala(perms(iperm_2+1,4)))
                                                elseif (iperm_2 == 5) then
                                                    W2 = &
CONJG(ktea2*wsq1(rq4_2,lala(perms(iperm_2+1,1)),lala(perms(iperm_2+1,2)),lala(perms(iperm_2+1,3)),lala(perms(iperm_2+1,4))))
                                                elseif (iperm_2 == 6) then
                                                    W2 = &
CONJG(ktea2*wsq1(rq4_2,lala(perms(iperm_2+1,1)),lala(perms(iperm_2+1,2)),lala(perms(iperm_2+1,3)),lala(perms(iperm_2+1,4))))
                                                elseif (iperm_2 == 7) then
                                                    W2 = &
ktea2*wsq1(rq4_2,lala(perms(iperm_2+1,1)),lala(perms(iperm_2+1,2)),lala(perms(iperm_2+1,3)),lala(perms(iperm_2+1,4)))
                                                end if
                                                wsqn(rq4_0,mu1,mu2,mu3,mu4) = &
wsqn(rq4_0,mu1,mu2,mu3,mu4) + 0.5d0*F(q5_1,q6_1,mu5,mu6)*W1*W2
                                            end if
                                            end do
                                            end if
                                            end do
                                            end if
                                            end do
                                            end do
                                            end if
                                            end do
                                            end do
                                        end if
                                        end do
                                        end do
                                        end if
                                        end do
                                        end do
                                        end if
                                        end do
                                        end do
                                        end if
                                        end do
                                        end do
                                        !$omp end do
                                        !$omp end parallel
                                    end if
                                end if
                            end do  
                        end do
                    end if
                end do
            end do
        end do
        wsq_scf = wsq_scf*(1-alpha_mix) +alpha_mix*(wsq1+wsqn)
        if (verbose) then
		    print "(A, ES16.6)", &
"     MAXVAL(|W(n)-W(n-1)|)             =", MAXVAL(ABS(wsq_scf-wsq2))
            print "(A, ES16.6)", &
"     Threshold (weighted by alpha_mix) =", eps*alpha_mix
            print*, ""
        end if
        if (MAXVAL(ABS(wsq_scf-wsq2)) < eps*alpha_mix) exit do1
        if (iter == maxiter) then
            wsq_scf = wsq1
            its_conv = .false.
        end if
    end do do1
    if (verbose) then
        tend = omp_get_wtime()
        if (its_conv) then
            print "(A, I4, A)", " Convergence found with", iter, " iterations"
        else
            print "(A, I4, A)", " Convergence not found with", iter, " iterations"
            print "(A)", " Returning first order W^{(1)}(q1,q2,q3,q4)"
        end if
        print*, ""
        print "(A, ES16.6,A)", "  Elapsed time inside get_scf_wsq():", tend-tstart, " seconds"
		print*, "=======================     DONE     ======================="
        print*, ""
	endif

end subroutine get_scf_wsq

subroutine get_sum_wsq(wsq1, F, refq4, refq4o, norbitq4, P, degs, verbose, &
wsq_scf, nrefq4, n_mode, iq, dimq4, nsym)
    
    implicit none

    complex(8), dimension(nrefq4,n_mode,n_mode,n_mode,n_mode), intent(in) :: wsq1

    double precision, dimension(iq, iq, n_mode, n_mode), intent(in) :: F

	integer, dimension(nrefq4,dimq4,4), intent(in) :: refq4
    integer, dimension(nrefq4,dimq4,2), intent(in) :: refq4o
    integer, dimension(nrefq4), intent(in) :: norbitq4

	complex(8), dimension(iq, n_mode, n_mode, nsym), intent(in) :: P
    logical, dimension(iq, n_mode, n_mode), intent(in) :: degs
    logical, intent(in) :: verbose

    complex(8), dimension(nrefq4,n_mode,n_mode,n_mode,n_mode), intent(out) :: wsq_scf

    complex(8), dimension(nrefq4,n_mode,n_mode,n_mode,n_mode) :: wsq2
	complex(8), dimension(nrefq4,n_mode,n_mode,n_mode,n_mode) :: wsqn

    double precision :: tstart, tend

    integer, dimension(12,4) :: perms
    integer, dimension(4) :: nunu, lala, mumu1, mumu2
    integer :: nrefq4, n_mode, iq, dimq4, nsym
    integer :: mu1,mu2,mu3,mu4,mu5,mu6,iter
    integer :: rq4_0,rq4_1,rq4_2,q4i_1,q4i_2,q1,q2,q3,q4,q5_1,q6_1,q5_2,q6_2
	integer :: q1_ref, q2_ref, q3_ref, q4_ref, q5_1ref, q5_2ref, q6_1ref, q6_2ref
	integer :: nu1,nu2,nu5,nu6,lambda3,lambda4,lambda5,lambda6
    integer :: q10,q20,q30,q40,q5,q6,iperm_1,iperm_2,isym_1,isym_2

	complex(8) :: ktea1, ktea2, W1, W2
    double precision, parameter :: eps = 1e-6
    integer, parameter :: maxiter = 50

    perms = reshape([ &
        1,2,3,4, 2,1,4,3, 3,4,1,2, 4,3,2,1, &
        1,3,2,4, 3,1,4,2, 2,4,1,3, 4,2,3,1, &
        4,2,3,1, 2,4,1,3, 3,1,4,2, 1,3,2,4], &
    [12,4], order=[2,1])

    if (verbose) then
        tstart = omp_get_wtime()
        print*, "======================= get_sum_wsq() ======================="
        print*, "Looking where to truncate the analytical sum of \Theta(-q1,q2,q3,-q4)."
        print*, ""
    end if
	wsq_scf = wsq1
    do1 : do iter = 1, maxiter
        if (verbose) then
            print*, "    Iteration", iter
            print*, "    -----------------------"
        end if
        if (iter == 1) then 
            wsq2 = wsq1
        else
            wsq2 = wsqn
        end if
        wsqn = (0.0d0,0.0d0)
        do rq4_0 = 1, nrefq4
            q10 = refq4(rq4_0,1,1)+1
            q20 = refq4(rq4_0,1,2)+1
            q30 = refq4(rq4_0,1,3)+1
            q40 = refq4(rq4_0,1,4)+1
            do rq4_1 = 1, nrefq4
                q1_ref = refq4(rq4_1,1,1)+1
                q2_ref = refq4(rq4_1,1,2)+1
                q5_1ref = refq4(rq4_1,1,3)+1
                q6_1ref = refq4(rq4_1,1,4)+1
                do q4i_1 = 1, norbitq4(rq4_1)
                    q1 = refq4(rq4_1,q4i_1,1)+1
                    q2 = refq4(rq4_1,q4i_1,2)+1
                    q5_1 = refq4(rq4_1,q4i_1,3)+1
                    q6_1 = refq4(rq4_1,q4i_1,4)+1
                    iperm_1 = refq4o(rq4_1,q4i_1,1)
					isym_1 = refq4o(rq4_1,q4i_1,2)+1
                    if ((q10 == q1) .and. (q20 == q2)) then
                        do rq4_2 = 1, nrefq4	
                            q5_2ref = refq4(rq4_2,1,1)+1
                            q6_2ref = refq4(rq4_2,1,2)+1
                            q3_ref = refq4(rq4_2,1,3)+1
                            q4_ref = refq4(rq4_2,1,4)+1	
                            do q4i_2 = 1, norbitq4(rq4_2)
                                q5_2 = refq4(rq4_2,q4i_2,1)+1
                                q6_2 = refq4(rq4_2,q4i_2,2)+1
                                q3 = refq4(rq4_2,q4i_2,3)+1
                                q4 = refq4(rq4_2,q4i_2,4)+1
                                iperm_2 = refq4o(rq4_2,q4i_2,1)
								isym_2 = refq4o(rq4_2,q4i_2,2)+1
                                if ((q30 == q3) .and. (q40 == q4)) then
                                    if ((q5_1 == q5_2) .and. (q6_1 == q6_2)) then
                                        !$omp parallel &
                                        !$omp private (mu2,mu3,mu4,mu5,mu6) &
                                        !$omp private (nu1,nu2,nu5,nu6,lambda5,lambda6,lambda3,lambda4,nunu,lala,mumu1,mumu2)
                                        !$omp do schedule (dynamic, 1) private (ktea1,ktea2,W1,W2)         
                                        do mu1 = 1, n_mode
                                        do nu1 = 1, n_mode
                                        if (degs(q1,nu1,mu1)) then
                                        do mu2 = 1, n_mode
                                        do nu2 = 1, n_mode
                                        if (degs(q2,nu2,mu2)) then
                                        do mu5 = 1, n_mode
                                        do nu5 = 1, n_mode
                                        if (degs(q5_1,nu5,mu5)) then
                                        do mu6 = 1, n_mode
                                        mumu1 = [mu1,mu2,mu5,mu6]
                                        do nu6 = 1, n_mode
                                        nunu = [nu1,nu2,nu5,nu6]
                                        if (degs(q6_1,nu6,mu6)) then
                                            ktea1 = &
P(q1_ref,mumu1(perms(iperm_1+1,1)),nunu(perms(iperm_1+1,1)),isym_1)*&
CONJG(P(q2_ref,mumu1(perms(iperm_1+1,2)),nunu(perms(iperm_1+1,2)),isym_1))*&
CONJG(P(q5_1ref,mumu1(perms(iperm_1+1,3)),nunu(perms(iperm_1+1,3)),isym_1))*&
P(q6_1ref,mumu1(perms(iperm_1+1,4)),nunu(perms(iperm_1+1,4)),isym_1)
                                            if (iperm_1 == 0) then
                                                W1 = &
ktea1*wsq2(rq4_1,nunu(perms(iperm_1+1,1)),nunu(perms(iperm_1+1,2)),nunu(perms(iperm_1+1,3)),nunu(perms(iperm_1+1,4)))
                                            elseif (iperm_1 == 1) then
                                                W1 = &
CONJG(ktea1*wsq2(rq4_1,nunu(perms(iperm_1+1,1)),nunu(perms(iperm_1+1,2)),nunu(perms(iperm_1+1,3)),nunu(perms(iperm_1+1,4))))
                                            elseif (iperm_1 == 2) then
                                                W1 = &
CONJG(ktea1*wsq2(rq4_1,nunu(perms(iperm_1+1,1)),nunu(perms(iperm_1+1,2)),nunu(perms(iperm_1+1,3)),nunu(perms(iperm_1+1,4))))
                                            elseif (iperm_1 == 3) then
                                                W1 = &
ktea1*wsq2(rq4_1,nunu(perms(iperm_1+1,1)),nunu(perms(iperm_1+1,2)),nunu(perms(iperm_1+1,3)),nunu(perms(iperm_1+1,4)))
                                            elseif (iperm_1 == 4) then
                                                W1 = &
ktea1*wsq2(rq4_1,nunu(perms(iperm_1+1,1)),nunu(perms(iperm_1+1,2)),nunu(perms(iperm_1+1,3)),nunu(perms(iperm_1+1,4)))
                                            elseif (iperm_1 == 5) then
                                                W1 = &
CONJG(ktea1*wsq2(rq4_1,nunu(perms(iperm_1+1,1)),nunu(perms(iperm_1+1,2)),nunu(perms(iperm_1+1,3)),nunu(perms(iperm_1+1,4))))
                                            elseif (iperm_1 == 6) then
                                                W1 = &
CONJG(ktea1*wsq2(rq4_1,nunu(perms(iperm_1+1,1)),nunu(perms(iperm_1+1,2)),nunu(perms(iperm_1+1,3)),nunu(perms(iperm_1+1,4))))
                                            elseif (iperm_1 == 7) then
                                                W1 = &
ktea1*wsq2(rq4_1,nunu(perms(iperm_1+1,1)),nunu(perms(iperm_1+1,2)),nunu(perms(iperm_1+1,3)),nunu(perms(iperm_1+1,4)))
                                            elseif (iperm_1 == 8) then
                                                W1 = &
ktea1*wsq2(rq4_1,nunu(perms(iperm_1+1,1)),nunu(perms(iperm_1+1,2)),nunu(perms(iperm_1+1,3)),nunu(perms(iperm_1+1,4)))
                                            elseif (iperm_1 == 9) then
                                                W1 = &
CONJG(ktea1*wsq2(rq4_1,nunu(perms(iperm_1+1,1)),nunu(perms(iperm_1+1,2)),nunu(perms(iperm_1+1,3)),nunu(perms(iperm_1+1,4))))
                                            elseif (iperm_1 == 10) then
                                                W1 = &
CONJG(ktea1*wsq2(rq4_1,nunu(perms(iperm_1+1,1)),nunu(perms(iperm_1+1,2)),nunu(perms(iperm_1+1,3)),nunu(perms(iperm_1+1,4))))
                                            elseif (iperm_1 == 11) then
                                                W1 = &
ktea1*wsq2(rq4_1,nunu(perms(iperm_1+1,1)),nunu(perms(iperm_1+1,2)),nunu(perms(iperm_1+1,3)),nunu(perms(iperm_1+1,4)))
                                            end if
                                            do mu3 = 1, n_mode
                                            do lambda3 = 1, n_mode
                                            if (degs(q3,lambda3,mu3)) then
                                            do mu4 = 1, n_mode
                                            do lambda4 = 1, n_mode
                                            mumu2 = [mu5,mu6,mu3,mu4]
                                            if (degs(q4,lambda4,mu4)) then
                                            do lambda5 = 1, n_mode
                                            if (degs(q5_2,lambda5,mu5)) then
                                            do lambda6 = 1, n_mode
                                            lala = [lambda5,lambda6,lambda3,lambda4]
                                            if (degs(q6_2,lambda6,mu6)) then
                                                ktea2 = &
P(q5_2ref,mumu2(perms(iperm_2+1,1)),lala(perms(iperm_2+1,1)),isym_2)*&
CONJG(P(q6_2ref,mumu2(perms(iperm_2+1,2)),lala(perms(iperm_2+1,2)),isym_2))*&
CONJG(P(q3_ref,mumu2(perms(iperm_2+1,3)),lala(perms(iperm_2+1,3)),isym_2))*&
P(q4_ref,mumu2(perms(iperm_2+1,4)),lala(perms(iperm_2+1,4)),isym_2)
                                                if (iperm_2 == 0) then
                                                    W2 = &
ktea2*wsq1(rq4_2,lala(perms(iperm_2+1,1)),lala(perms(iperm_2+1,2)),lala(perms(iperm_2+1,3)),lala(perms(iperm_2+1,4)))
                                                elseif (iperm_2 == 1) then
                                                    W2 = &
CONJG(ktea2*wsq1(rq4_2,lala(perms(iperm_2+1,1)),lala(perms(iperm_2+1,2)),lala(perms(iperm_2+1,3)),lala(perms(iperm_2+1,4))))
                                                elseif (iperm_2 == 2) then
                                                    W2 = &
CONJG(ktea2*wsq1(rq4_2,lala(perms(iperm_2+1,1)),lala(perms(iperm_2+1,2)),lala(perms(iperm_2+1,3)),lala(perms(iperm_2+1,4))))
                                                elseif (iperm_2 == 3) then
                                                    W2 = &
ktea2*wsq1(rq4_2,lala(perms(iperm_2+1,1)),lala(perms(iperm_2+1,2)),lala(perms(iperm_2+1,3)),lala(perms(iperm_2+1,4)))
                                                elseif (iperm_2 == 4) then
                                                    W2 = &
ktea2*wsq1(rq4_2,lala(perms(iperm_2+1,1)),lala(perms(iperm_2+1,2)),lala(perms(iperm_2+1,3)),lala(perms(iperm_2+1,4)))
                                                elseif (iperm_2 == 5) then
                                                    W2 = &
CONJG(ktea2*wsq1(rq4_2,lala(perms(iperm_2+1,1)),lala(perms(iperm_2+1,2)),lala(perms(iperm_2+1,3)),lala(perms(iperm_2+1,4))))
                                                elseif (iperm_2 == 6) then
                                                    W2 = &
CONJG(ktea2*wsq1(rq4_2,lala(perms(iperm_2+1,1)),lala(perms(iperm_2+1,2)),lala(perms(iperm_2+1,3)),lala(perms(iperm_2+1,4))))
                                                elseif (iperm_2 == 7) then
                                                    W2 = &
ktea2*wsq1(rq4_2,lala(perms(iperm_2+1,1)),lala(perms(iperm_2+1,2)),lala(perms(iperm_2+1,3)),lala(perms(iperm_2+1,4)))
                                                elseif (iperm_2 == 8) then
                                                    W2 = &
ktea2*wsq1(rq4_2,lala(perms(iperm_2+1,1)),lala(perms(iperm_2+1,2)),lala(perms(iperm_2+1,3)),lala(perms(iperm_2+1,4)))
                                                elseif (iperm_2 == 9) then
                                                    W2 = &
CONJG(ktea2*wsq1(rq4_2,lala(perms(iperm_2+1,1)),lala(perms(iperm_2+1,2)),lala(perms(iperm_2+1,3)),lala(perms(iperm_2+1,4))))
                                                elseif (iperm_2 == 10) then
                                                    W2 = &
CONJG(ktea2*wsq1(rq4_2,lala(perms(iperm_2+1,1)),lala(perms(iperm_2+1,2)),lala(perms(iperm_2+1,3)),lala(perms(iperm_2+1,4))))
                                                elseif (iperm_2 == 11) then
                                                    W2 = &
ktea2*wsq1(rq4_2,lala(perms(iperm_2+1,1)),lala(perms(iperm_2+1,2)),lala(perms(iperm_2+1,3)),lala(perms(iperm_2+1,4)))
                                                end if
                                                wsqn(rq4_0,mu1,mu2,mu3,mu4) = &
wsqn(rq4_0,mu1,mu2,mu3,mu4) + 0.5d0*F(q5_1,q6_1,mu5,mu6)*W1*W2
                                            end if
                                            end do
                                            end if
                                            end do
                                            end if
                                            end do
                                            end do
                                            end if
                                            end do
                                            end do
                                        end if
                                        end do
                                        end do
                                        end if
                                        end do
                                        end do
                                        end if
                                        end do
                                        end do
                                        end if
                                        end do
                                        end do
                                        !$omp end do
                                        !$omp end parallel
                                    end if
                                end if
                            end do  
                        end do
                    end if
                end do
            end do
        end do
        if (verbose) then
		    print*, "    MAXVAL(|W^{n}|)=", MAXVAL(ABS(wsqn)), ". Threshold=", eps
            print*, ""
        end if
		wsq_scf = wsq_scf + wsqn
        if (MAXVAL(ABS(wsqn)) < eps) exit do1
        if (MAXVAL(ABS(wsqn)) > MAXVAL(ABS(wsq2))) then
            print*, "W^{(n)} is not decaying. Returning first order correction..."
            wsqn=wsq1
        end if
    end do do1
    if (verbose) then
        tend = omp_get_wtime()
        print*, "Convergence found with", iter, "iterations"
        print*, ""
        print "(A, ES16.6,A)", "  Elapsed time inside get_scf_wsq():", tend-tstart, " seconds"
		print*, "=======================     DONE     ======================="
        print*, ""
	endif

end subroutine get_sum_wsq

subroutine get_indep2fc( &
        vs_red, refq2, refq2o, norbitq2, orbit2t, n_indep_elem, indep_elem, rot_cart, mapping, map_uc, map_tr, T_list, q_list, F, &
        verbose, indep_fc, nrefq2, n_mode, dimq2, nref2, dim2, Natom_sc, Nqpoint, Nsym)

	implicit none

    complex, dimension(nrefq2,n_mode,n_mode,n_mode) :: vs_red
    integer, dimension(nrefq2,dimq2,2), intent(in) :: refq2
    integer, dimension(nrefq2,dimq2,2), intent(in) :: refq2o
    integer, dimension(nrefq2), intent(in) :: norbitq2

    integer, dimension(nref2,dim2,2), intent(in) :: orbit2t
    integer, dimension(nref2), intent(in) :: n_indep_elem
    integer, dimension(nref2,9), intent(in) :: indep_elem

    double precision, dimension(Nsym,3,3), intent(in) :: rot_cart

    integer, dimension(Natom_sc, Nsym), intent(in) :: mapping
    integer, dimension(Natom_sc, Natom_sc), intent(in) :: map_uc
    integer, dimension(Natom_sc), intent(in) :: map_tr
    double precision, dimension(Nqpoint, 3), intent(in) :: T_list, q_list
    double precision, dimension(Nqpoint, Nqpoint, n_mode, n_mode), intent(in) :: F

    logical, intent(in) :: verbose

    double precision, dimension(nref2,9), intent(out) :: indep_fc

    double precision, parameter :: pi  = 4.0d0*atan(1.0d0)
    complex, parameter :: j = (0.0d0, 1.0d0)

    integer, dimension(2) :: munu
    integer, dimension(2,2) :: permutations
    integer :: ref2, nat1, nat2, i, alpha, beta, rq2, q2i, iperm, isym, nat1p, nat2p, nat1p_uc, nat2p_uc, mu, nu
    integer :: nrefq2, n_mode, dimq2, nref2, dim2, Natom_sc, Nqpoint, Nsym
    complex, dimension(n_mode) :: vs_q
    complex :: Vsa, Vsb

    if (verbose) then
        print*, "======================= get_indep2fc() ======================="
        print*, "Computing independent third order corrections to FCs \Phi^{3}_{a,b}^{\alpha,\beta}"
    end if

    permutations = reshape([1,2,2,1],[2,2])

    indep_fc = 0
    do ref2 = 1, nref2
        nat1 = orbit2t(ref2,1,1)
        nat2 = orbit2t(ref2,1,2)
        do i = 1, n_indep_elem(ref2)
            alpha = indep_elem(ref2,i)/3
            beta = mod(indep_elem(ref2,i),3)
            do rq2 = 1, nrefq2
                do q2i = 1, norbitq2(rq2)
                    iperm = refq2o(rq2,q2i,1)
                    isym = refq2o(rq2,q2i,2)
                    nat1p = findloc(mapping(:,isym+1),nat1, dim=1)-1
                    nat2p = findloc(mapping(:,isym+1),nat2, dim=1)-1
                    nat1p_uc = map_uc(nat1p+1,nat1p+1)
                    nat2p_uc = map_uc(nat2p+1,nat2p+1)
                    do mu = 1, n_mode
                        do nu = 1, n_mode
                            munu(1) = mu
                            munu(2) = nu
                            vs_q = &
vs_red(rq2,:,munu(permutations(iperm+1,1)), munu(permutations(iperm+1,2)))
                            Vsa = &
dot_product(rot_cart(isym+1,alpha+1,:),vs_q(nat1p_uc*3+1:(nat1p_uc+1)*3+1)) * &
exp(-2*j*pi*dot_product(q_list(refq2(rq2,1,1)+1,:)-q_list(refq2(rq2,1,2)+1,:),T_list(map_tr(nat1p+1)+1,:)))
                            Vsb = &
dot_product(rot_cart(isym+1,beta+1,:),vs_q(nat2p_uc*3+1:(nat2p_uc+1)*3+1)) * &
exp(-2*j*pi*dot_product(q_list(refq2(rq2,1,1)+1,:)-q_list(refq2(rq2,1,2)+1,:),T_list(map_tr(nat2p+1)+1,:)))
                            if (iperm == 0) then
                                indep_fc(ref2,i) = indep_fc(ref2,i) + &
real(0.5d0*F(refq2(rq2,1,1)+1,refq2(rq2,1,2)+1,munu(permutations(iperm+1,1)), munu(permutations(iperm+1,2)))*Vsa*CONJG(Vsb))
                            else 
                                indep_fc(ref2,i) = indep_fc(ref2,i) + &
real(0.5d0*F(refq2(rq2,1,1)+1,refq2(rq2,1,2)+1,munu(permutations(iperm+1,1)), munu(permutations(iperm+1,2)))*CONJG(Vsa)*Vsb)
                            end if
                        end do
                    end do
                end do
            end do
        end do
    end do
    if (verbose) then
        print *, "=======================     DONE     ======================="
        print *, ""
    end if
end subroutine get_indep2fc

subroutine get_indep2fc_v4( &
    vs, ws_red, refq4, refq4o, norbitq4, orbit2t, n_indep_elem, indep_elem, & 
    F, Ds, degs, mapping, rot_cart, verbose, indep_fc, nref2, nrefq4, n_mode, dimq4, dim2, n_mode_sc, iq, nsym, nat_sc)

	implicit none

    complex, dimension(n_mode_sc,iq,iq,n_mode,n_mode), intent(in) :: vs
    complex, dimension(nrefq4,n_mode,n_mode,n_mode,n_mode), intent(in) :: ws_red
    integer, dimension(nrefq4,dimq4,4), intent(in) :: refq4
    integer, dimension(nrefq4,dimq4,2), intent(in) :: refq4o
    integer, dimension(nrefq4), intent(in) :: norbitq4

    integer, dimension(nref2,dim2,2), intent(in) :: orbit2t
    integer, dimension(nref2), intent(in) :: n_indep_elem
    integer, dimension(nref2,9), intent(in) :: indep_elem
    double precision, dimension(iq, iq, n_mode, n_mode), intent(in) :: F

    complex(8), dimension(iq, n_mode, n_mode, nsym), intent(in) :: Ds
    logical, dimension(iq, n_mode, n_mode), intent(in) :: degs
    
    integer, dimension(nat_sc, nsym), intent(in) :: mapping
    double precision, dimension(nsym,3,3), intent(in) :: rot_cart

    logical, intent(in) :: verbose

    double precision, dimension(nref2,9), intent(out) :: indep_fc

    double precision, parameter :: pi  = 4.0d0*atan(1.0d0)
    complex, parameter :: j = (0.0d0, 1.0d0)

    integer, dimension(4) :: munu,nunu
    integer :: ref2, nat1, nat2, i, alpha, beta, rq4, q4i, iperm, isym, nsym
    integer :: mu1, mu2, mu3, mu4, nu1, nu2, nu3, nu4, q1, q2, q3, q4, q10, q20, q30, q40
    integer :: nrefq4, n_mode, dimq4, nref2, dim2, n_mode_sc, iq, nat1p, nat2p, nat_sc
    complex :: W, vsa, vsb

    if (verbose) then
        print*, "======================= get_indep2fc_v4() ======================="
        print*, "Computing independent 4th order corrections to FCs \Phi^{4}_{a,b}^{\alpha,\beta}"
    end if

    indep_fc = 0.0d0
    !$omp parallel &
    !$omp private (nat1,nat2,i,alpha,beta,rq4,q4i,q1,q2,q3,q4,q10,q20,q30,q40) &
    !$omp private (iperm,isym,nat1p,nat2p,mu1,mu2,mu3,mu4)
    !$omp do schedule (dynamic, 1) private (W,vsa,vsb) 
    do ref2 = 1, nref2
        nat1 = orbit2t(ref2,1,1)
        nat2 = orbit2t(ref2,1,2)
        do i = 1, n_indep_elem(ref2)
            alpha = indep_elem(ref2,i)/3+1
            beta = mod(indep_elem(ref2,i),3)+1
            do rq4 = 1, nrefq4
                do q4i = 1, norbitq4(rq4)
                    q1 = refq4(rq4,q4i,1)+1
                    q2 = refq4(rq4,q4i,2)+1
                    q3 = refq4(rq4,q4i,3)+1
                    q4 = refq4(rq4,q4i,4)+1
                    q10 = refq4(rq4,1,1)+1
                    q20 = refq4(rq4,1,2)+1
                    q30 = refq4(rq4,1,3)+1
                    q40 = refq4(rq4,1,4)+1
                    iperm = refq4o(rq4,q4i,1)
                    isym = refq4o(rq4,q4i,2)
                    nat1p = findloc(mapping(:,isym+1),nat1, dim=1)-1
                    nat2p = findloc(mapping(:,isym+1),nat2, dim=1)-1
                    do mu1=1,n_mode
                    do mu2=1,n_mode
                    do mu3=1,n_mode
                    do mu4=1,n_mode
                        W = ws_red(rq4,mu1,mu2,mu3,mu4)
                        if (iperm == 0) then
                            vsa= &
    dot_product(rot_cart(isym+1,alpha,:),vs(3*nat1p+1:3*(nat1p+1)+1,q10,q20,mu1,mu2))
                            vsb= &
    dot_product(rot_cart(isym+1,beta,:),vs(3*nat2p+1:3*(nat2p+1)+1,q30,q40,mu3,mu4))
                            indep_fc(ref2,i) = indep_fc(ref2,i) + real(0.25d0* &
F(q10,q20,mu1,mu2) * &
F(q30,q40,mu3,mu4) * &
vsa*W*CONJG(vsb))
                        elseif (iperm == 1) then
                            vsa= &
    dot_product(rot_cart(isym+1,alpha,:),vs(3*nat1p+1:3*(nat1p+1)+1,q20,q10,mu2,mu1))
                            vsb= &
    dot_product(rot_cart(isym+1,beta,:),vs(3*nat2p+1:3*(nat2p+1)+1,q40,q30,mu4,mu3))
                            indep_fc(ref2,i) = indep_fc(ref2,i) + real(0.25d0* &
F(q10,q20,mu1,mu2) * &
F(q30,q40,mu3,mu4) * &
vsa*CONJG(W)*CONJG(vsb))
                        elseif (iperm == 2) then
                            vsa= &
    dot_product(rot_cart(isym+1,alpha,:),vs(3*nat1p+1:3*(nat1p+1)+1,q30,q40,mu3,mu4))
                            vsb= &
    dot_product(rot_cart(isym+1,beta,:),vs(3*nat2p+1:3*(nat2p+1)+1,q10,q20,mu1,mu2))
                            indep_fc(ref2,i) = indep_fc(ref2,i) + real(0.25d0* &
F(q10,q20,mu1,mu2) * &
F(q30,q40,mu3,mu4) * &
vsa*CONJG(W)*CONJG(vsb))
                        elseif (iperm == 3) then
                            vsa= &
    dot_product(rot_cart(isym+1,alpha,:),vs(3*nat1p+1:3*(nat1p+1)+1,q40,q30,mu4,mu3))
                            vsb= &
    dot_product(rot_cart(isym+1,beta,:),vs(3*nat2p+1:3*(nat2p+1)+1,q20,q10,mu2,mu1))
                            indep_fc(ref2,i) = indep_fc(ref2,i) + real(0.25d0* &
F(q10,q20,mu1,mu2) * &
F(q30,q40,mu3,mu4) * &
vsa*W*CONJG(vsb))
                        elseif (iperm == 4) then
                            vsa= &
    dot_product(rot_cart(isym+1,alpha,:),vs(3*nat1p+1:3*(nat1p+1)+1,q10,q30,mu1,mu3))
                            vsb= &
    dot_product(rot_cart(isym+1,beta,:),vs(3*nat2p+1:3*(nat2p+1)+1,q20,q40,mu2,mu4))
                            indep_fc(ref2,i) = indep_fc(ref2,i) + real(0.25d0* &
F(q10,q20,mu1,mu2) * &
F(q30,q40,mu3,mu4) * &
vsa*W*CONJG(vsb))
                        elseif (iperm == 5) then
                            vsa= &
    dot_product(rot_cart(isym+1,alpha,:),vs(3*nat1p+1:3*(nat1p+1)+1,q30,q10,mu3,mu1))
                            vsb= &
    dot_product(rot_cart(isym+1,beta,:),vs(3*nat2p+1:3*(nat2p+1)+1,q40,q20,mu4,mu2))
                            indep_fc(ref2,i) = indep_fc(ref2,i) + real(0.25d0* &
F(q10,q20,mu1,mu2) * &
F(q30,q40,mu3,mu4) * &
vsa*CONJG(W)*CONJG(vsb))
                        elseif (iperm == 6) then
                            vsa= &
    dot_product(rot_cart(isym+1,alpha,:),vs(3*nat1p+1:3*(nat1p+1)+1,q20,q40,mu2,mu4))
                            vsb= &
    dot_product(rot_cart(isym+1,beta,:),vs(3*nat2p+1:3*(nat2p+1)+1,q10,q30,mu1,mu3))
                            indep_fc(ref2,i) = indep_fc(ref2,i) + real(0.25d0* &
F(q10,q20,mu1,mu2) * &
F(q30,q40,mu3,mu4) * &
vsa*CONJG(W)*CONJG(vsb))
                        elseif (iperm == 7) then
                            vsa= &
    dot_product(rot_cart(isym+1,alpha,:),vs(3*nat1p+1:3*(nat1p+1)+1,q40,q20,mu4,mu2))
                            vsb= &
    dot_product(rot_cart(isym+1,beta,:),vs(3*nat2p+1:3*(nat2p+1)+1,q30,q10,mu3,mu1))
                            indep_fc(ref2,i) = indep_fc(ref2,i) + real(0.25d0* &
F(q10,q20,mu1,mu2) * &
F(q30,q40,mu3,mu4) * &
vsa*W*CONJG(vsb))
                        end if
                    end do
                    end do
                    end do
                    end do
                end do 
            end do
        end do
    end do
    !$omp end do
    !$omp end parallel
    if (verbose) then
        print *, "=======================     DONE     ======================="
        print *, ""
    end if
end subroutine get_indep2fc_v4


subroutine get_rotq_ws(q, pol_vecs, v4, T_list, q_list, map_tr, map_uc, isym, rot_cart_inv, rot_4fc_inv, mapping, &
    w, iq, n_mode, n_mode_sc, nat_sc, nsym)

    implicit none

    integer, dimension(4), intent(in) :: q
    complex(8), dimension(iq,n_mode,n_mode_sc), intent(in) :: pol_vecs
    double precision, dimension(n_mode_sc,n_mode_sc,n_mode_sc,n_mode_sc), intent(in) :: v4
    double precision, dimension(iq, 3), intent(in) :: T_list, q_list
    
    integer, dimension(nat_sc), intent(in) :: map_tr
    integer, dimension(nat_sc, nat_sc), intent(in) :: map_uc

    integer, intent(in) :: isym
    double precision, dimension(nsym,3,3), intent(in) :: rot_cart_inv
    double precision, dimension(24,nsym,81,81), intent(in) :: rot_4fc_inv
    integer, dimension(iq,nsym), intent(in) :: mapping
    !integer, dimension(nat_sc,nsym), intent(in) :: mapping_inv

    complex(8), dimension(n_mode,n_mode,n_mode,n_mode), intent(out) :: w
    double precision, dimension(81) :: flatv4, rot_flat_v4

    complex(8) :: rot_pol_alpha1, rot_pol_alpha2, rot_pol_alpha3, rot_pol_alpha4
    integer :: iq, n_mode, n_mode_sc, nat_sc

    integer :: q1,q2,q3,q4,mu1,mu2,mu3,mu4,flatmu,index,ii
    integer :: a1,a2,a3,a4,alpha1,alpha2,alpha3,alpha4,nsym
    
    double precision, parameter :: pi  = 4.0d0*atan(1.0d0)
    complex, parameter :: j = (0.0d0, 1.0d0)
    complex :: phase

    ! Get integers
    w = 0.0d0

    q1 = q(1)+1
    q2 = q(2)+1
    q3 = q(3)+1
    q4 = q(4)+1

    !$omp parallel private (mu1,mu2,mu3,mu4, rot_pol_alpha1, rot_pol_alpha2, rot_pol_alpha3, rot_pol_alpha4)
    !$omp do schedule (dynamic, 1) private (a1,a2,a3,a4,alpha1,alpha2,alpha3,alpha4,ii,rot_flat_v4,flatv4,index)
    do flatmu = 1, n_mode**4
        mu1 = (flatmu-1)/n_mode**3 + 1
        mu2 = mod(flatmu-1,n_mode**3)/n_mode**2 + 1
        mu3 = mod(flatmu-1,n_mode**2)/n_mode + 1
        mu4 = mod(flatmu-1,n_mode) + 1
        do a1 = 1, nat_sc
            do a2 = 1, nat_sc
                do a3 = 1, nat_sc
                    do a4 = 1, nat_sc
                        do alpha1 = 1, 3
                            do alpha2 = 1, 3
                                do alpha3 = 1, 3
                                    do alpha4 = 1, 3
                                        index = 3*(3*(3*(alpha1-1)+alpha2-1)+alpha3-1)+alpha4
                                        flatv4(index) = &
v4(3*(a1-1)+alpha1,3*(a2-1)+alpha2,3*(a3-1)+alpha3,3*(a4-1)+alpha4)
                                    end do
                                end do
                            end do
                        end do
                        do alpha1 = 1, 3
                            rot_pol_alpha1 = &
! dot_product(rot_cart_inv(isym+1,alpha1,:),pol_vecs(mapping(q1,isym+1)+1,mu1,3*(a1-1)+1:3*a1+1))
dot_product(rot_cart_inv(1,alpha1,:),pol_vecs(mapping(q1,isym+1)+1,mu1,3*(a1-1)+1:3*a1+1))
                            do alpha2 = 1, 3
                                rot_pol_alpha2 = &
! dot_product(rot_cart_inv(isym+1,alpha2,:),pol_vecs(mapping(q2,isym+1)+1,mu2,3*(a2-1)+1:3*a2+1))
dot_product(rot_cart_inv(1,alpha2,:),pol_vecs(mapping(q2,isym+1)+1,mu2,3*(a2-1)+1:3*a2+1))
                                do alpha3 = 1, 3
                                    rot_pol_alpha3 = &
! dot_product(rot_cart_inv(isym+1,alpha3,:),pol_vecs(mapping(q3,isym+1)+1,mu3,3*(a3-1)+1:3*a3+1))
dot_product(rot_cart_inv(1,alpha3,:),pol_vecs(mapping(q3,isym+1)+1,mu3,3*(a3-1)+1:3*a3+1))
                                    do alpha4 = 1, 3
                                        rot_pol_alpha4 = &
! dot_product(rot_cart_inv(isym+1,alpha4,:),pol_vecs(mapping(q4,isym+1)+1,mu4,3*(a4-1)+1:3*a4+1))
dot_product(rot_cart_inv(1,alpha4,:),pol_vecs(mapping(q4,isym+1)+1,mu4,3*(a4-1)+1:3*a4+1))
                                        index = 3*(3*(3*(alpha1-1)+alpha2-1)+alpha3-1)+alpha4
                                        rot_flat_v4(index) = 0
                                        do ii = 1, 81
                                            rot_flat_v4(index) = rot_flat_v4(index) + &
        ! rot_4fc_inv(1, isym+1, index, ii) * flatv4(ii)
        rot_4fc_inv(1, 1, index, ii) * flatv4(ii)
                                        end do
                                        w(mu1,mu2,mu3,mu4) = w(mu1,mu2,mu3,mu4) + &
conjg(rot_pol_alpha1*rot_pol_alpha2)*&
rot_flat_v4(index)*&
rot_pol_alpha3*rot_pol_alpha4
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
    !$omp end do
    !$omp end parallel
    print*, q1, q2, q3, q4, MAXVAL(ABS(w))

end subroutine get_rotq_ws

  ! This subroutine calculates the stochastic average of an element of the
! third order force constants based on fu2 and lmat. 
! The force constants element is given with three indexes that correspond 
! to atom-cartesian coordinates.
!
! NOTE: there is a first addend that is not implemented at the moment
!       because it should be zero if the odd correction is calculated
!       at the atomic positions that minimize the gradient with respect
!       to Wyckoff positions. This non-implemented addend vanishes
!       if positions are fixed by symmetry.
! Version de get_v3_element_sym pero guardando datos en una matrix auxiliar

subroutine get_v3_element_sym (na_in, nb_in, nc_in, ur, eprod, f, rho, log_err, &
        nsym, s_cart, irt, translations_irt, v3, nat_sc, n_mode, n_random, nr)

    use omp_lib
    use stochastic
    implicit none

    integer, intent(in) :: na_in, nb_in, nc_in, nsym
    double precision, dimension(n_random,n_mode), intent(in) :: ur
    double precision, dimension(n_mode,n_mode), intent(in) :: eprod
    double precision, dimension(n_random,nat_sc,3), intent(in) :: f
    double precision, dimension(n_random), intent(in) :: rho
    character (len=10), intent(in) :: log_err
    double precision, dimension(3,3,48), intent(in) :: s_cart
    integer, dimension(48,nat_sc), intent(in) :: irt
    integer, dimension(nat_sc,nr), intent(in) :: translations_irt
    double precision, intent(out) :: v3

    integer :: nat_sc, n_mode, n_random, nr
  
    double precision, dimension(:), allocatable :: fun 
    double precision, dimension(:,:), allocatable :: Rot,v3_tp, v3_s!aux matrix where the elements have transl and perm sym applied
    double precision :: v3_aux, av, av_err
    double precision :: aux_elem1, aux_elem2, aux_elem3
    
    integer, dimension(3) :: triplet, triplet_sym
    integer :: alpha, beta, gamma 
    integer :: isym, ii, jj, kk, na, nb, nc, cartesian_index, r
  
    logical, parameter :: debug = .false.
    real :: t1, t2

    ! Allocate stuff
    if (debug) then
      print *, "=== V3 ELEMENT DEBUG ==="
      print *, "Atomic-cartesian coord of the element:", na_in, nb_in, nc_in
      print *, ""
      call flush()
    end if
  
    allocate(fun(n_random))
    allocate(Rot(nsym,27))
    allocate(v3_tp(nsym,27))
    allocate(v3_s(nsym,nsym))

    v3 = 0.0d0
    
    triplet(1) = na_in/3 + 1 ! +1 Fortran
    triplet(2) = nb_in/3 + 1
    triplet(3) = nc_in/3 + 1 
    alpha = mod(na_in,3) + 1 
    beta = mod(nb_in,3) + 1
    gamma = mod(nc_in,3) + 1

    v3_tp(:,:) = 0.0d0
    do isym = 1, nsym
        triplet_sym(1) = irt(isym, triplet(1))
        triplet_sym(2) = irt(isym, triplet(2))
        triplet_sym(3) = irt(isym, triplet(3))
        do ii = 1, 3
            do jj = 1, 3
                do kk = 1, 3
                        ! Single atomic-cartesian index in fortran
                        na = 3 * (triplet_sym(1)-1) + ii
                        nb = 3 * (triplet_sym(2)-1) + jj
                        nc = 3 * (triplet_sym(3)-1) + kk
                        cartesian_index = (3*(ii-1)+(jj-1))*3+(kk-1)+1 

                        ! Obtain Rot element with trans and perm symmetries applied only if the multiplication
                        ! of Rotation matrix elements is not 0:
                        Rot(isym,cartesian_index) = s_cart(alpha, ii, isym)*s_cart(beta, jj, isym)*s_cart(gamma, kk, isym)
                        if (Rot(isym,cartesian_index) /= 0.0d0 .and. abs(Rot(isym,cartesian_index))>1.0e-12) then

                                !v3_aux = 0.0d0
                                if (na == nb.and.na == nc.and.nb == nc) then !No need of applying perm symmetry to the element:
                                        do r = 1, nr
                                               ! Translated Single atomic-cartesian index (nb, nc same as na)
                                               na = 3 * (translations_irt(triplet_sym(1), r)-1) + ii
                                               fun(:) = ( eprod(na,na) - ur(:,na) * ur(:,na) ) * &
                                                       f(:, translations_irt(triplet_sym(1), r), ii)
                                               call average_error_weight(fun,rho,log_err,av,av_err)
                                               v3_tp(isym, cartesian_index) = v3_tp(isym, cartesian_index) + av
                                        end do

                                elseif (na == nb .and. na /= nc) then
                                        ! Need to apply perm+trans symmetry to the element
                                        do r = 1, nr
                                               ! Translated Single atomic-cartesian index in fortran
                                               na = 3 * (translations_irt(triplet_sym(1), r)-1) + ii
                                               nc = 3 * (translations_irt(triplet_sym(3), r)-1) + kk
                                               fun(:) = ( eprod(na,nc) - ur(:,na) * ur(:,nc) ) * &
                                                       f(:, translations_irt(triplet_sym(1), r), ii)
                                               call average_error_weight(fun,rho,log_err,aux_elem1,av_err)
                                               fun(:) = ( eprod(na,na) - ur(:,na) * ur(:,na) ) * &
                                                       f(:, translations_irt(triplet_sym(3), r), kk)
                                               call average_error_weight(fun,rho,log_err,aux_elem2,av_err)
                                               v3_tp(isym, cartesian_index) = v3_tp(isym, cartesian_index) + &
                                                       (2.0d0*aux_elem1+aux_elem2) / 3.0d0
                                        end do
                                                                
                                elseif (na == nc .and. na /= nb) then
                                        do r = 1, nr
                                               ! Translated Single atomic-cartesian index in fortran
                                               na = 3 * (translations_irt(triplet_sym(1), r)-1) + ii
                                               nb = 3 * (translations_irt(triplet_sym(2), r)-1) + jj
                                               fun(:) = ( eprod(nb,na) - ur(:,nb) * ur(:,na) ) * &
                                                       f(:, translations_irt(triplet_sym(1), r), ii)
                                               call average_error_weight(fun,rho,log_err,aux_elem1,av_err)
                                               fun(:) = ( eprod(na,na) - ur(:,na) * ur(:,na) ) * &
                                                       f(:, translations_irt(triplet_sym(2), r), jj)
                                               call average_error_weight(fun,rho,log_err,aux_elem2,av_err)
                                               v3_tp(isym, cartesian_index) = v3_tp(isym, cartesian_index) + &
                                                       (2.0d0*aux_elem1+aux_elem2) / 3.0d0
                                        end do

                                elseif (nb == nc .and. na /= nc) then
                                        do r = 1, nr
                                               ! Translated Single atomic-cartesian index in fortran
                                               na = 3 * (translations_irt(triplet_sym(1), r)-1) + ii
                                               nb = 3 * (translations_irt(triplet_sym(2), r)-1) + jj
                                               fun(:) = ( eprod(nb,nb) - ur(:,nb) * ur(:,nb) ) * &
                                                       f(:, translations_irt(triplet_sym(1), r), ii)
                                               call average_error_weight(fun,rho,log_err,aux_elem1,av_err)
                                               fun(:) = ( eprod(na,nb) - ur(:,na) * ur(:,nb) ) * &
                                                       f(:, translations_irt(triplet_sym(2), r), jj)
                                               call average_error_weight(fun,rho,log_err,aux_elem2,av_err)
                                               v3_tp(isym, cartesian_index) = v3_tp(isym, cartesian_index) + &
                                                       (aux_elem1+2.0d0*aux_elem2) / 3.0d0
                                        end do

                                elseif (na /= nb .and. na /= nc .and. nb /= nc) then
                                        do r = 1, nr
                                               ! Translated Single atomic-cartesian index in fortran
                                               na = 3 * (translations_irt(triplet_sym(1), r)-1) + ii
                                               nb = 3 * (translations_irt(triplet_sym(2), r)-1) + jj
                                               nc = 3 * (translations_irt(triplet_sym(3), r)-1) + kk
                                               ! fun is inv with respect nb <-> nc
                                               fun(:) = ( eprod(nb,nc) - ur(:,nb) * ur(:,nc) ) * &
                                                       f(:, translations_irt(triplet_sym(1), r), ii)
                                               call average_error_weight(fun,rho,log_err,aux_elem1,av_err)
                                               fun(:) = ( eprod(na,nc) - ur(:,na) * ur(:,nc) ) * &
                                                       f(:, translations_irt(triplet_sym(2), r), jj)
                                               call average_error_weight(fun,rho,log_err,aux_elem2,av_err)
                                               fun(:) = ( eprod(na,nb) - ur(:,na) * ur(:,nb) ) * &
                                                       f(:, translations_irt(triplet_sym(3), r), kk)
                                               call average_error_weight(fun,rho,log_err,aux_elem3,av_err)
                                               v3_tp(isym, cartesian_index) = v3_tp(isym, cartesian_index) + &
                                                       (aux_elem1+aux_elem2+aux_elem3) / 3.0d0
                                        end do
                                end if
                        end if        
                end do
            end do
        end do
    end do
   
    ! Apply space group symmetry     
    !call dgemm('N', 'T', nsym, nsym, 27, 1.0d0, &
    !        v3_tp, nsym, Rot, 27, 0.0d0, v3_s, nsym)
  
    do isym = 1, nsym
        v3 = v3 + sum(Rot(isym,:)*v3_tp(isym,:))
    end do
    v3 = v3 / (real(nr)*real(nsym))

    deallocate(fun,v3_tp,Rot,v3_s)
  
  end subroutine get_v3_element_sym

subroutine get_v4_element_sym (na_in, nb_in, nc_in, nd_in, ur, eprod, f, rho, log_err, &
        nsym, s_cart, irt, translations_irt, v4, nat_sc, n_mode, n_random, nr)

        use omp_lib
        use stochastic
        implicit none

        integer, intent(in) :: na_in, nb_in, nc_in, nd_in, nsym
        double precision, dimension(n_random,n_mode), intent(in) :: ur
        double precision, dimension(n_mode,n_mode), intent(in) :: eprod
        double precision, dimension(n_random,nat_sc,3), intent(in) :: f
        double precision, dimension(n_random), intent(in) :: rho
        character (len=10), intent(in) :: log_err
        double precision, dimension(3,3,48), intent(in) :: s_cart
        integer, dimension(48,nat_sc), intent(in) :: irt
        integer, dimension(nat_sc,nr), intent(in) :: translations_irt
        double precision, intent(out) :: v4

        integer :: nat_sc, n_mode, n_random, nr

        double precision, dimension(:), allocatable :: fun 
        double precision, dimension(:,:), allocatable :: Rot,v4_tp, v4_s!aux matrix where the elements have transl and perm sym applied
        double precision :: v4_aux, av, av_err
        double precision :: aux_elem1, aux_elem2, aux_elem3, aux_elem4, aux_elem5, aux_elem6

        integer, dimension(4) :: qplet, qplet_sym
        integer :: alpha, beta, gamma, delta
        integer :: isym, ii, jj, kk, ll, na, nb, nc, nd, cartesian_index, r

        logical, parameter :: debug = .false.
        real :: t1, t2

        ! Allocate stuff
        if (debug) then
        print *, "=== V4 ELEMENT DEBUG ==="
        print *, "Atomic-cartesian coord of the element:", na_in, nb_in, nc_in, nd_in
        print *, ""
        call flush()
        end if

        allocate(fun(n_random))
        allocate(Rot(nsym,81))
        allocate(v4_tp(nsym,81))
        allocate(v4_s(nsym,nsym))

        v4 = 0.0d0

        qplet(1) = na_in/3 + 1
        qplet(2) = nb_in/3 + 1
        qplet(3) = nc_in/3 + 1 
        qplet(4) = nd_in/3 + 1 
        alpha = mod(na_in,3) + 1 
        beta = mod(nb_in,3) + 1
        gamma = mod(nc_in,3) + 1
        delta = mod(nd_in,3) + 1

        v4_tp(:,:) = 0.0d0
        do isym = 1, nsym
                qplet_sym(1) = irt(isym, qplet(1))
                qplet_sym(2) = irt(isym, qplet(2))
                qplet_sym(3) = irt(isym, qplet(3))
                qplet_sym(4) = irt(isym, qplet(4))
                do ii = 1, 3
                        do jj = 1, 3
                                do kk = 1, 3
                                        do ll = 1, 3
                                                ! Single atomic-cartesian index in fortran
                                                na = 3 * (qplet_sym(1)-1) + ii
                                                nb = 3 * (qplet_sym(2)-1) + jj
                                                nc = 3 * (qplet_sym(3)-1) + kk
                                                nd = 3 * (qplet_sym(4)-1) + ll
                                                cartesian_index = (((ii-1)*3+(jj-1))*3+(kk-1))*3+ll
                                                Rot(isym,cartesian_index) = &
                s_cart(alpha, ii, isym)*s_cart(beta, jj, isym)*s_cart(gamma, kk, isym)*s_cart(delta, ll, isym)

                                                if (abs(Rot(isym,cartesian_index))>1.0e-12) then !If element is zero, we need nothing.
                                                        if (na == nb .and. na == nc .and. na == nd) then ! All elements are the same.
                                                                !print*, "1"
                                                                do r = 1, nr
                                                                        na = 3 * (translations_irt(qplet_sym(1), r)-1) + ii
                                                                        fun(:) = &
                - ur(:,na) * ur(:,na) * ur(:,na) * f(:,translations_irt(qplet_sym(1), r), ii)
                                                                        call average_error_weight(fun,rho,log_err,av,av_err)
                                                                        v4_tp(isym, cartesian_index) = &
                v4_tp(isym, cartesian_index) + av
                                                                end do
                                                        elseif (na == nb .and. na == nc) then ! 3 elements are equal.
                                                                !print*, "2.1"
                                                                do r = 1, nr
                                                                        na = 3 * (translations_irt(qplet_sym(1), r)-1) + ii
                                                                        nd = 3 * (translations_irt(qplet_sym(4), r)-1) + ll
                                                                        fun(:) = &
                - ur(:,na) * ur(:,na) * ur(:,na) * f(:,translations_irt(qplet_sym(4), r), ll)
                                                                        call average_error_weight(fun,rho,log_err,aux_elem1,av_err)
                                                                        fun(:) = &
                - ur(:,nd) * ur(:,na) * ur(:,na) * f(:,translations_irt(qplet_sym(1), r), ii)
                                                                        call average_error_weight(fun,rho,log_err,aux_elem2,av_err)
                                                                        v4_tp(isym, cartesian_index) = &
                v4_tp(isym, cartesian_index) + (aux_elem1+3*aux_elem2)/4.0d0
                                                                end do
                                                        elseif (nb == nc .and. nb == nd) then ! 3 elements are equal.
                                                                !print*, "2.2"
                                                                do r = 1, nr
                                                                        nb = 3 * (translations_irt(qplet_sym(2), r)-1) + jj
                                                                        na = 3 * (translations_irt(qplet_sym(1), r)-1) + ii
                                                                        fun(:) = &
                - ur(:,nb) * ur(:,nb) * ur(:,nb) * f(:,translations_irt(qplet_sym(1), r), ii)
                                                                        call average_error_weight(fun,rho,log_err,aux_elem1,av_err)
                                                                        fun(:) = &
                - ur(:,na) * ur(:,nb) * ur(:,nb) * f(:,translations_irt(qplet_sym(2), r), jj)
                                                                        call average_error_weight(fun,rho,log_err,aux_elem2,av_err)
                                                                        v4_tp(isym, cartesian_index) = &
                v4_tp(isym, cartesian_index) + (aux_elem1+3*aux_elem2)/4.0d0
                                                                end do
                                                        elseif (nc == nd .and. nc == na) then ! 3 elements are equal.
                                                                !print*, "2.3"
                                                                do r = 1, nr
                                                                        nc = 3 * (translations_irt(qplet_sym(3), r)-1) + kk
                                                                        nb = 3 * (translations_irt(qplet_sym(2), r)-1) + jj
                                                                        fun(:) = &
                - ur(:,nc) * ur(:,nc) * ur(:,nc) * f(:,translations_irt(qplet_sym(2), r), jj)
                                                                        call average_error_weight(fun,rho,log_err,aux_elem1,av_err)
                                                                        fun(:) = &
                - ur(:,nb) * ur(:,nc) * ur(:,nc) * f(:,translations_irt(qplet_sym(3), r), kk)
                                                                        call average_error_weight(fun,rho,log_err,aux_elem2,av_err)
                                                                        v4_tp(isym, cartesian_index) = &
                v4_tp(isym, cartesian_index) + (aux_elem1+3*aux_elem2)/4.0d0
                                                                end do
                                                        elseif (nd == na .and. nd == nb) then ! 3 elements are equal.
                                                                !print*, "2.4"
                                                                do r = 1, nr
                                                                        nd = 3 * (translations_irt(qplet_sym(4), r)-1) + ll
                                                                        nc = 3 * (translations_irt(qplet_sym(3), r)-1) + kk
                                                                        fun(:) = &
                - ur(:,nd) * ur(:,nd) * ur(:,nd) * f(:,translations_irt(qplet_sym(3), r), kk)
                                                                        call average_error_weight(fun,rho,log_err,aux_elem1,av_err)
                                                                        fun(:) = &
                - ur(:,nc) * ur(:,nd) * ur(:,nd) * f(:,translations_irt(qplet_sym(4), r), ll)
                                                                        call average_error_weight(fun,rho,log_err,aux_elem2,av_err)
                                                                        v4_tp(isym, cartesian_index) = &
                v4_tp(isym, cartesian_index) + (aux_elem1+3*aux_elem2)/4.0d0
                                                                end do
                                                        elseif (na == nb .and. nc == nd) then ! 2 and 2
                                                                !print*, "3.1"
                                                                do r = 1, nr
                                                                        na = 3 * (translations_irt(qplet_sym(1), r)-1) + ii
                                                                        nc = 3 * (translations_irt(qplet_sym(3), r)-1) + kk
                                                                        fun(:) = &
                - ur(:,na) * ur(:,na) * ur(:,nc) * f(:,translations_irt(qplet_sym(3), r), kk)                                               
                                                                        call average_error_weight(fun,rho,log_err,aux_elem1,av_err)
                                                                        fun(:) = &
                - ur(:,nc) * ur(:,nc) * ur(:,na) * f(:,translations_irt(qplet_sym(1), r), ii)                                               
                                                                        call average_error_weight(fun,rho,log_err,aux_elem2,av_err)
                                                                        v4_tp(isym, cartesian_index) = &
                v4_tp(isym, cartesian_index) + (aux_elem1+aux_elem2) / 2.0d0
                                                                end do
                                                        elseif ((na == nc .and. nb == nd) .or. (na==nd .and. nb==nc)) then
                                                                !print*, "3.2"
                                                                do r = 1, nr
                                                                        na = 3 * (translations_irt(qplet_sym(1), r)-1) + ii
                                                                        nb = 3 * (translations_irt(qplet_sym(2), r)-1) + jj
                                                                        fun(:) = &
                - ur(:,na) * ur(:,na) * ur(:,nb) * f(:,translations_irt(qplet_sym(2), r), jj)                                               
                                                                        call average_error_weight(fun,rho,log_err,aux_elem1,av_err)
                                                                        fun(:) = &
                - ur(:,nb) * ur(:,nb) * ur(:,na) * f(:,translations_irt(qplet_sym(1), r), ii)
                                                                        call average_error_weight(fun,rho,log_err,aux_elem2,av_err)
                                                                        v4_tp(isym, cartesian_index) = &
                v4_tp(isym, cartesian_index) + (aux_elem1+aux_elem2) / 2.0d0
                                                                end do
                                                        elseif (na == nb) then ! 2-1
                                                                !print*, "4.1"
                                                                do r = 1, nr
                                                                        na = 3 * (translations_irt(qplet_sym(1), r)-1) + ii
                                                                        nc = 3 * (translations_irt(qplet_sym(3), r)-1) + kk
                                                                        nd = 3 * (translations_irt(qplet_sym(4), r)-1) + ll
                                                                        fun(:) = &
                - ur(:,na) * ur(:,na) * ur(:,nc) * f(:,translations_irt(qplet_sym(4), r), ll)                                               
                                                                        call average_error_weight(fun,rho,log_err,aux_elem1,av_err)
                                                                        fun(:) = &
                - ur(:,na) * ur(:,na) * ur(:,nd) * f(:,translations_irt(qplet_sym(3), r), kk)                                               
                                                                        call average_error_weight(fun,rho,log_err,aux_elem2,av_err)
                                                                        fun(:) = &
                - ur(:,na) * ur(:,nc) * ur(:,nd) * f(:,translations_irt(qplet_sym(1), r), ii)                                               
                                                                        call average_error_weight(fun,rho,log_err,aux_elem3,av_err)
                                                                        v4_tp(isym, cartesian_index) = &
                v4_tp(isym, cartesian_index) + (aux_elem1+aux_elem2+2*aux_elem3)/4.0d0
                                                                end do
                                                        elseif (na == nc) then ! 2-1
                                                                !print*, "4.2"
                                                                do r = 1, nr
                                                                        na = 3 * (translations_irt(qplet_sym(1), r)-1) + ii
                                                                        nb = 3 * (translations_irt(qplet_sym(2), r)-1) + jj
                                                                        nd = 3 * (translations_irt(qplet_sym(4), r)-1) + ll
                                                                        fun(:) = &
                - ur(:,na) * ur(:,nb) * ur(:,na) * f(:,translations_irt(qplet_sym(4), r), ll)                                               
                                                                        call average_error_weight(fun,rho,log_err,aux_elem1,av_err)
                                                                        fun(:) = &
                - ur(:,na) * ur(:,na) * ur(:,nd) * f(:,translations_irt(qplet_sym(2), r), jj)                                               
                                                                        call average_error_weight(fun,rho,log_err,aux_elem2,av_err)
                                                                        fun(:) = &
                - ur(:,na) * ur(:,nb) * ur(:,nd) * f(:,translations_irt(qplet_sym(1), r), ii)                                               
                                                                        call average_error_weight(fun,rho,log_err,aux_elem3,av_err)
                                                                        v4_tp(isym, cartesian_index) = &
                v4_tp(isym, cartesian_index) + (aux_elem1+aux_elem2+2*aux_elem3)/4.0d0
                                                                end do
                                                        elseif (na == nd) then ! 2-1
                                                                !print*, "4.3"
                                                                do r = 1, nr
                                                                        na = 3 * (translations_irt(qplet_sym(1), r)-1) + ii
                                                                        nb = 3 * (translations_irt(qplet_sym(2), r)-1) + jj
                                                                        nc = 3 * (translations_irt(qplet_sym(3), r)-1) + kk
                                                                        fun(:) = &
                - ur(:,na) * ur(:,nb) * ur(:,na) * f(:,translations_irt(qplet_sym(3), r), kk)                                               
                                                                        call average_error_weight(fun,rho,log_err,aux_elem1,av_err)
                                                                        fun(:) = &
                - ur(:,na) * ur(:,na) * ur(:,nc) * f(:,translations_irt(qplet_sym(2), r), jj)                                               
                                                                        call average_error_weight(fun,rho,log_err,aux_elem2,av_err)
                                                                        fun(:) = &
                - ur(:,na) * ur(:,nb) * ur(:,nc) * f(:,translations_irt(qplet_sym(1), r), ii)                                               
                                                                        call average_error_weight(fun,rho,log_err,aux_elem3,av_err)
                                                                        v4_tp(isym, cartesian_index) = &
                v4_tp(isym, cartesian_index) + (aux_elem1+aux_elem2+2*aux_elem3)/4.0d0
                                                                end do
                                                        elseif (nb == nc) then ! 2-1
                                                                !print*, "4.4"
                                                                do r = 1, nr
                                                                        na = 3 * (translations_irt(qplet_sym(1), r)-1) + ii
                                                                        nb = 3 * (translations_irt(qplet_sym(2), r)-1) + jj
                                                                        nd = 3 * (translations_irt(qplet_sym(4), r)-1) + ll
                                                                        fun(:) = &
                - ur(:,na) * ur(:,nb) * ur(:,nb) * f(:,translations_irt(qplet_sym(4), r), ll)                                               
                                                                        call average_error_weight(fun,rho,log_err,aux_elem1,av_err)
                                                                        fun(:) = &
                - ur(:,nb) * ur(:,nb) * ur(:,nd) * f(:,translations_irt(qplet_sym(1), r), ii)                                               
                                                                        call average_error_weight(fun,rho,log_err,aux_elem2,av_err)
                                                                        fun(:) = &
                - ur(:,na) * ur(:,nb) * ur(:,nd) * f(:,translations_irt(qplet_sym(2), r), jj)                                               
                                                                        call average_error_weight(fun,rho,log_err,aux_elem3,av_err)
                                                                        v4_tp(isym, cartesian_index) = &
                v4_tp(isym, cartesian_index) + (aux_elem1+aux_elem2+2*aux_elem3)/4.0d0
                                                                end do
                                                        elseif (nb == nd) then ! 2-1
                                                                !print*, "4.5"
                                                                do r = 1, nr
                                                                        na = 3 * (translations_irt(qplet_sym(1), r)-1) + ii
                                                                        nb = 3 * (translations_irt(qplet_sym(2), r)-1) + jj
                                                                        nc = 3 * (translations_irt(qplet_sym(3), r)-1) + kk
                                                                        fun(:) = &
                - ur(:,na) * ur(:,nb) * ur(:,nb) * f(:,translations_irt(qplet_sym(3), r), kk)                                               
                                                                        call average_error_weight(fun,rho,log_err,aux_elem1,av_err)
                                                                        fun(:) = &
                - ur(:,nb) * ur(:,nb) * ur(:,nc) * f(:,translations_irt(qplet_sym(1), r), ii)                                               
                                                                        call average_error_weight(fun,rho,log_err,aux_elem2,av_err)
                                                                        fun(:) = &
                - ur(:,na) * ur(:,nb) * ur(:,nc) * f(:,translations_irt(qplet_sym(2), r), jj)                                               
                                                                        call average_error_weight(fun,rho,log_err,aux_elem3,av_err)
                                                                        v4_tp(isym, cartesian_index) = &
                v4_tp(isym, cartesian_index) + (aux_elem1+aux_elem2+2*aux_elem3)/4.0d0
                                                                end do
                                                        elseif (nc == nd) then ! 2-1
                                                                !print*, "4.6"
                                                                do r = 1, nr
                                                                        na = 3 * (translations_irt(qplet_sym(1), r)-1) + ii
                                                                        nb = 3 * (translations_irt(qplet_sym(2), r)-1) + jj
                                                                        nc = 3 * (translations_irt(qplet_sym(3), r)-1) + kk
                                                                        fun(:) = &
                - ur(:,na) * ur(:,nc) * ur(:,nc) * f(:,translations_irt(qplet_sym(2), r), jj)                                               
                                                                        call average_error_weight(fun,rho,log_err,aux_elem1,av_err)
                                                                        fun(:) = &
                - ur(:,nb) * ur(:,nc) * ur(:,nc) * f(:,translations_irt(qplet_sym(1), r), ii)                                               
                                                                        call average_error_weight(fun,rho,log_err,aux_elem2,av_err)
                                                                        fun(:) = &
                - ur(:,na) * ur(:,nb) * ur(:,nc) * f(:,translations_irt(qplet_sym(3), r), kk)                                               
                                                                        call average_error_weight(fun,rho,log_err,aux_elem3,av_err)
                                                                        v4_tp(isym, cartesian_index) = &
                v4_tp(isym, cartesian_index) + (aux_elem1+aux_elem2+2*aux_elem3) / 4.0d0
                                                                        end do               
                                                        else ! None
                                                                do r = 1, nr
                                                                        na = 3 * (translations_irt(qplet_sym(1), r)-1) + ii
                                                                        nb = 3 * (translations_irt(qplet_sym(2), r)-1) + jj
                                                                        nc = 3 * (translations_irt(qplet_sym(3), r)-1) + kk
                                                                        nd = 3 * (translations_irt(qplet_sym(4), r)-1) + ll
                                                                        fun(:) = &
                - ur(:,na) * ur(:,nb) * ur(:,nc) * f(:,translations_irt(qplet_sym(4), r), ll)                                               
                                                                        call average_error_weight(fun,rho,log_err,aux_elem1,av_err)
                                                                        fun(:) = &
                - ur(:,nb) * ur(:,nc) * ur(:,nd) * f(:,translations_irt(qplet_sym(1), r), ii)                                               
                                                                        call average_error_weight(fun,rho,log_err,aux_elem2,av_err)
                                                                        fun(:) = &
                - ur(:,na) * ur(:,nc) * ur(:,nd) * f(:,translations_irt(qplet_sym(2), r), jj)                                               
                                                                        call average_error_weight(fun,rho,log_err,aux_elem3,av_err)
                                                                        fun(:) = &
                - ur(:,na) * ur(:,nb) * ur(:,nd) * f(:,translations_irt(qplet_sym(3), r), kk)                                               
                                                                        call average_error_weight(fun,rho,log_err,aux_elem4,av_err)
                                                                        v4_tp(isym, cartesian_index) = &
                v4_tp(isym, cartesian_index) + (aux_elem1+aux_elem2+aux_elem3+aux_elem4) / 4.0d0
                                                                end do
                                                        end if
                                                end if        
                                        end do
                                end do
                        end do
                end do
        end do

        do isym = 1, nsym
                v4 = v4 + sum(Rot(isym,:)*v4_tp(isym,:))
        end do
        v4 = v4 / (real(nr)*real(nsym))

        deallocate(fun,v4_tp,Rot,v4_s)
  
end subroutine get_v4_element_sym

subroutine get_nref3(nat, nat_sc, tot3, nsym, mappings, map_uc, nref3)

        integer, intent(in) :: nat, nat_sc, nsym, tot3
        integer, dimension(:,:), intent(in) :: mappings, map_uc

        integer, intent(out) :: nref3

        integer, dimension(tot3,3) :: all3, ref3
        integer, dimension(6*nsym, 3) :: equilist

        integer, dimension(6,3) :: permutations
        integer, dimension(3) :: triplet, triplet_perm, triplet_sym
        integer :: ii, jj, kk, nall3, equiv, iperm, isym
        logical :: its_in_list

        permutations = reshape( &
        [1,2,3, 2,1,3, 3,2,1, 1,3,2, 2,3,1, 3,1,2], &
        [6,3], order=[2,1])

        nref3 = 0
        nall3 = 0
        doii : do ii = 0, nat-1
                dojj : do jj = 0, nat_sc-1
                        dokk: do kk = 0, nat_sc-1
                                triplet = [ii, jj, kk]
                                call triplet_in_list(triplet, all3, nall3, its_in_list)
                                if (its_in_list) cycle dokk
                                nref3 = nref3 + 1
                                ref3(nref3,:) = triplet
                                equiv = 0
                                do iperm = 1, 6
                                        triplet_perm(1) = triplet(permutations(iperm,1))
                                        triplet_perm(2) = triplet(permutations(iperm,2))
                                        triplet_perm(3) = triplet(permutations(iperm,3))
                                        do isym = 1, nsym
                                                triplet_sym = &
[mappings(triplet_perm(1)+1,isym), mappings(triplet_perm(2)+1,isym), mappings(triplet_perm(3)+1,isym)]
                                                if (triplet_sym(1) >= nat) then
                                                        triplet_sym(2) = map_uc(triplet_sym(1)+1, triplet_sym(2)+1)
                                                        triplet_sym(3) = map_uc(triplet_sym(1)+1, triplet_sym(3)+1)
                                                        triplet_sym(1) = map_uc(triplet_sym(1)+1, triplet_sym(1)+1)
                                                end if
                                                call triplet_in_list(triplet_sym, equilist, equiv, its_in_list)
                                                if (((iperm==1) .and. (isym==1)) .or. (.not. its_in_list)) then
                                                        equilist(equiv+1,:) = triplet_sym
                                                        all3(nall3+1,:) = triplet_sym
                                                        equiv = equiv + 1
                                                        nall3 = nall3 + 1
                                                end if
                                        end do
                                end do
                                if (nall3 == tot3) exit doii
                        end do dokk
                end do dojj
        end do doii
end subroutine

subroutine get_nref4(nat, nat_sc, tot4, nsym, mappings, map_uc, nref4)

        integer, intent(in) :: nat, nat_sc, nsym, tot4
        integer, dimension(:,:), intent(in) :: mappings, map_uc

        integer, intent(out) :: nref4

        integer, dimension(tot4,4) :: all4, ref4
        integer, dimension(24*nsym, 4) :: equilist

        integer, dimension(24,4) :: permutations
        integer, dimension(4) :: qplet, qplet_perm, qplet_sym
        integer :: ii, jj, kk, ll, nall4, equiv, iperm, isym
        logical :: its_in_list

        permutations = reshape( &
        [1,2,3,4, 2,1,3,4, 3,2,1,4, 1,3,2,4, 2,3,1,4, 3,1,2,4, &
        1,2,4,3, 2,1,4,3, 4,2,1,3, 1,4,2,3, 2,4,1,3, 4,1,2,3, &
        1,4,3,2, 4,1,3,2, 3,4,1,2, 1,3,4,2, 4,3,1,2, 3,1,4,2, &
        4,2,3,1, 2,4,3,1, 3,2,4,1, 4,3,2,1, 2,3,4,1, 3,4,2,1], &
        [24,4], order=[2,1])

        nref4 = 0
        nall4 = 0
        doii : do ii = 0, nat-1
                dojj : do jj = 0, nat_sc-1
                        dokk: do kk = 0, nat_sc-1
                                doll: do ll = 0, nat_sc-1
                                        qplet = [ii, jj, kk, ll]
                                        call qplet_in_list(qplet, all4, nall4, its_in_list)
                                        if (its_in_list) cycle doll
                                        nref4 = nref4 + 1
                                        ref4(nref4,:) = qplet
                                        equiv = 0
                                        do iperm = 1, 24
                                                qplet_perm(1) = qplet(permutations(iperm,1))
                                                qplet_perm(2) = qplet(permutations(iperm,2))
                                                qplet_perm(3) = qplet(permutations(iperm,3))
                                                qplet_perm(4) = qplet(permutations(iperm,4))
                                                do isym = 1, nsym
                                                        qplet_sym = &
[mappings(qplet_perm(1)+1,isym), mappings(qplet_perm(2)+1,isym), mappings(qplet_perm(3)+1,isym), mappings(qplet_perm(4)+1,isym)]
                                                        if (qplet_sym(1) >= nat) then
                                                                qplet_sym(2) = map_uc(qplet_sym(1)+1, qplet_sym(2)+1)
                                                                qplet_sym(3) = map_uc(qplet_sym(1)+1, qplet_sym(3)+1)
                                                                qplet_sym(4) = map_uc(qplet_sym(1)+1, qplet_sym(4)+1)
                                                                qplet_sym(1) = map_uc(qplet_sym(1)+1, qplet_sym(1)+1)
                                                        end if
                                                        call qplet_in_list(qplet_sym, equilist, equiv, its_in_list)
                                                        if (((iperm==1) .and. (isym==1)) .or. (.not. its_in_list)) then
                                                                equilist(equiv+1,:) = qplet_sym
                                                                all4(nall4+1,:) = qplet_sym
                                                                equiv = equiv + 1
                                                                nall4 = nall4 + 1
                                                        end if
                                                end do
                                        end do
                                        if (nall4 == tot4) exit doii
                                end do doll
                        end do dokk
                end do dojj
        end do doii
end subroutine


subroutine get_q_nref4(mappings, nref4, norbitq4, iq, nsym)

        integer, dimension(iq,nsym), intent(in) :: mappings

        integer, intent(out) :: nref4
        integer, intent(out) :: norbitq4

        integer, allocatable, dimension(:,:) :: all4
        integer, allocatable, dimension(:) :: norbit4
        integer, dimension(12*nsym, 4) :: equilist

        integer, dimension(8,4) :: permutations
        integer, dimension(4) :: qplet, qplet_perm, qplet_sym
        integer :: q1, q2, q3, q4, nall4, equiv, iperm, isym, iq, nsym
        logical :: its_in_list

        allocate(all4(iq**4,4))
        allocate(norbit4(iq**4))

        permutations = reshape([ &
            1,2,3,4, 2,1,4,3, 3,4,1,2, 4,3,2,1, &
            1,3,2,4, 3,1,4,2, 2,4,1,3, 4,2,3,1], &
        [8,4], order=[2,1])

        !permutations = reshape([1,2,3,4], [1,4])

        nref4 = 0
        nall4 = 0
        norbit4 = 0
        do1 : do q1 = 1, iq
                do2 : do q2 = 1, iq
                        do3: do q3 = 1, iq
                                do4: do q4 = 1, iq
                                        qplet = [q1-1, q2-1, q3-1, q4-1]
                                        call qplet_in_list(qplet, all4, nall4, its_in_list)
                                        if (its_in_list) cycle do4
                                        nref4 = nref4 + 1
                                        equiv = 0
                                        do iperm = 1, 8
                                                qplet_perm(1) = qplet(permutations(iperm,1))
                                                qplet_perm(2) = qplet(permutations(iperm,2))
                                                qplet_perm(3) = qplet(permutations(iperm,3))
                                                qplet_perm(4) = qplet(permutations(iperm,4))
                                                do isym = 1, nsym
                                                        qplet_sym = &
[mappings(qplet_perm(1)+1,isym), mappings(qplet_perm(2)+1,isym), mappings(qplet_perm(3)+1,isym), mappings(qplet_perm(4)+1,isym)]
                                                        call qplet_in_list(qplet_sym, equilist, equiv, its_in_list)
                                                        if (((iperm==1) .and. (isym==1)) .or. (.not. its_in_list)) then
                                                                equilist(equiv+1,:) = qplet_sym
                                                                all4(nall4+1,:) = qplet_sym
                                                                equiv = equiv + 1
                                                                nall4 = nall4 + 1
                                                        end if
                                                end do
                                        end do
                                        norbit4(nref4) = equiv
                                        if (nall4 == iq**4) exit do1
                                end do do4
                        end do do3
                end do do2
        end do do1
        norbitq4 = MAXVAL(norbit4)
end subroutine

subroutine recognize_q_quadruplet(nref4, norbitq4, q_list, mappings, verbose, &
orbit4t, orbit4o, norbit, ref4, iq, nsym)
        integer, intent(in) :: nref4, norbitq4
        double precision, dimension(iq,3), intent(in) :: q_list
        integer, dimension(iq,nsym), intent(in) :: mappings
        logical, intent(in) :: verbose

        integer, dimension(nref4,norbitq4,4), intent(out) :: orbit4t
        integer, dimension(nref4,norbitq4,2), intent(out) :: orbit4o
        integer, dimension(nref4), intent(out) :: norbit
        integer, intent(out) :: ref4


        integer, dimension(iq**4,4) :: all4
        integer, dimension(norbitq4, 4) :: equilist

        integer, dimension(8,4) :: permutations
        integer, dimension(4) :: qplet, qplet_perm, qplet_sym
        integer :: q1, q2, q3, q4, nall4, equiv, iperm, isym
        integer :: iq, nsym

        double precision, dimension(3) :: qsum
        logical :: its_in_list, its_cons

        orbit4t = 0
        orbit4o = 0
        
        permutations = reshape([ &
            1,2,3,4, 2,1,4,3, 3,4,1,2, 4,3,2,1, &
            1,3,2,4, 3,1,4,2, 2,4,1,3, 4,2,3,1], &
        [8,4], order=[2,1])

        ref4 = 0
        nall4 = 0
        do1 : do q1 = 1, iq
                do2 : do q2 = 1, iq
                        do3: do q3 = 1, iq
                                do4: do q4 = 1, iq
                                        qplet = [q1-1,q2-1,q3-1,q4-1]
                                        call qplet_in_list(qplet, all4, nall4, its_in_list)
                                        if (its_in_list) cycle do4
                                        qsum(:) = &
-q_list(q1,:)+q_list(q2,:)+q_list(q3,:)-q_list(q4,:)
                                        call is_qcons(qsum, its_cons)
                                        if (its_cons) then
                                            ref4 = ref4 + 1                                        
                                            equiv = 0
                                            do iperm = 1, 8
                                                    qplet_perm(1) = qplet(permutations(iperm,1))
                                                    qplet_perm(2) = qplet(permutations(iperm,2))
                                                    qplet_perm(3) = qplet(permutations(iperm,3))
                                                    qplet_perm(4) = qplet(permutations(iperm,4))
                                                    do isym = 1, nsym
                                                            qplet_sym = &
    [mappings(qplet_perm(1)+1,isym), mappings(qplet_perm(2)+1,isym), mappings(qplet_perm(3)+1,isym), mappings(qplet_perm(4)+1,isym)]
                                                            call qplet_in_list(qplet_sym, equilist, equiv, its_in_list)
                                                            if (((iperm==1) .and. (isym==1)) .or. (.not. its_in_list)) then
                                                                    equilist(equiv+1,:) = qplet_sym
                                                                    orbit4t(ref4,equiv+1,:) = qplet_sym
                                                                    all4(nall4+1,:) = qplet_sym
                                                                    orbit4o(ref4,equiv+1,:) = [iperm-1, isym-1]
                                                                    equiv = equiv + 1
                                                                    nall4 = nall4 + 1
                                                            end if
                                                    end do
                                            end do
                                            norbit(ref4) = equiv
                                        else                                
                                            equiv = 0
                                            do iperm = 1, 8
                                                    qplet_perm(1) = qplet(permutations(iperm,1))
                                                    qplet_perm(2) = qplet(permutations(iperm,2))
                                                    qplet_perm(3) = qplet(permutations(iperm,3))
                                                    qplet_perm(4) = qplet(permutations(iperm,4))
                                                    do isym = 1, nsym
                                                            qplet_sym = &
    [mappings(qplet_perm(1)+1,isym), mappings(qplet_perm(2)+1,isym), mappings(qplet_perm(3)+1,isym), mappings(qplet_perm(4)+1,isym)]
                                                            call qplet_in_list(qplet_sym, equilist, equiv, its_in_list)
                                                            if (((iperm==1) .and. (isym==1)) .or. (.not. its_in_list)) then
                                                                    equilist(equiv+1,:) = qplet_sym
                                                                    all4(nall4+1,:) = qplet_sym
                                                                    equiv = equiv + 1
                                                                    nall4 = nall4 + 1
                                                            end if
                                                    end do
                                            end do
                                        end if
                                        if (verbose) then
                                            print*, "Reference q-quadruplet: (",q1,",",q2,",",q3,",",q4,")"
                                            print*, "Orbit size:", equiv
                                            print*, ""
                                        end if
                                        if (nall4 == iq**4) exit do1
                                end do do4
                        end do do3
                end do do2
        end do do1
end subroutine

subroutine recognize_triplet(nat, nat_sc, tot3, nref3, nsym, mappings, map_uc, nontrivial, M, &
    verbose, orbit3t, orbit3o, norbit, indep_fc, n_indep_fc, kernel, mapping_triplet)

        integer, intent(in) :: nat, nat_sc, nsym, tot3, nref3
        integer, dimension(:,:,:), intent(in) :: nontrivial
        double precision, dimension(:,:,:,:), intent(in) :: M
        integer, dimension(:,:), intent(in) :: mappings, map_uc
        logical, intent(in) :: verbose

        double precision, dimension(nref3,27,27), intent(out) :: kernel
        integer, dimension(nref3,27), intent(out) :: indep_fc 
        integer, dimension(nref3), intent(out) :: norbit, n_indep_fc
        integer, dimension(nref3,6*nsym,3), intent(out) :: orbit3t
        integer, dimension(nref3,6*nsym,2), intent(out) :: orbit3o
        integer, dimension(nat,nat_sc,nat_sc,3), intent(out) :: mapping_triplet

        double precision, allocatable, dimension(:,:) :: constrain_reduced

        integer, dimension(tot3,3) :: all3
        integer, dimension(6*nsym, 3) :: equilist
        integer, dimension(27) :: indep
        double precision, dimension(27,27) :: kern
        double precision, dimension(6*nsym*27,27) :: constrain

        integer, dimension(6,3) :: permutations
        integer, dimension(3) :: triplet, triplet_perm, triplet_sym
        integer :: ii, jj, kk, nall3, equiv, nconstrain, iperm, isym, iaux, &
jaux, indexprime, ll, nindep, iw, ref3
        logical :: its_in_list

        character(len=100) :: filename


        kernel = 0
        orbit3t = 0
        orbit3o = 0
        mapping_triplet = 0
        permutations = reshape( &
        [1,2,3, 2,1,3, 3,2,1, 1,3,2, 2,3,1, 3,1,2], &
        [6,3], order=[2,1])

        ref3 = 0
        nall3 = 0
        doii : do ii = 0, nat-1
                dojj : do jj = 0, nat_sc-1
                        dokk: do kk = 0, nat_sc-1
                                triplet = [ii, jj, kk]
                                call triplet_in_list(triplet, all3, nall3, its_in_list)
                                if (its_in_list) cycle dokk
                                ref3 = ref3 + 1
                                equiv = 0
                                nconstrain = 0
                                constrain = 0
                                do iperm = 1, 6
                                        triplet_perm(1) = triplet(permutations(iperm,1))
                                        triplet_perm(2) = triplet(permutations(iperm,2))
                                        triplet_perm(3) = triplet(permutations(iperm,3))
                                        do isym = 1, nsym
                                                triplet_sym = &
[mappings(triplet_perm(1)+1,isym), mappings(triplet_perm(2)+1,isym), mappings(triplet_perm(3)+1,isym)]
                                                if (triplet_sym(1) >= nat) then
                                                        triplet_sym(2) = map_uc(triplet_sym(1)+1, triplet_sym(2)+1)
                                                        triplet_sym(3) = map_uc(triplet_sym(1)+1, triplet_sym(3)+1)
                                                        triplet_sym(1) = map_uc(triplet_sym(1)+1, triplet_sym(1)+1)
                                                end if
                                                call triplet_in_list(triplet_sym, equilist, equiv, its_in_list)
                                                if (((iperm==1) .and. (isym==1)) .or. (.not. its_in_list)) then
                                                        equilist(equiv+1,:) = triplet_sym
                                                        orbit3t(ref3,equiv+1,:) = triplet_sym
                                                        all3(nall3+1,:) = triplet_sym
                                                        orbit3o(ref3,equiv+1,:) = [iperm-1, isym-1]
                                                        mapping_triplet(triplet_sym(1)+1,triplet_sym(2)+1,triplet_sym(3)+1,:) = &
[ref3-1,iperm-1,isym-1]
                                                        equiv = equiv + 1
                                                        nall3 = nall3 + 1
                                                end if
                                                if (all(triplet == triplet_sym)) then
                                                        do indexprime = 1, 27
                                                                if (nontrivial(iperm,isym,indexprime)==1) then
                                                                        do ll = 1, 27
                                                                                constrain(nconstrain+1, ll) = & 
                                                                        M(iperm,isym,indexprime,ll)
                                                                        end do
                                                                        nconstrain = nconstrain + 1
                                                                end if
                                                        end do
                                                end if
                                        end do
                                end do
                                norbit(ref3) = equiv
                                allocate(constrain_reduced(max(nconstrain,27),27))
                                constrain_reduced(:,:) = 0
                                do iaux = 1, nconstrain
                                        do jaux =1, 27
                                                constrain_reduced(iaux,jaux) = constrain(iaux,jaux)
                                        end do
                                end do
                                call gauss_jordan(constrain_reduced, max(nconstrain,27), 27, nconstrain, kern, indep, nindep)

                                deallocate(constrain_reduced)
                                if (verbose) then
                                    print *, "Reference triplet: (", ii, ",", jj, ",", kk, ")."
                                    print *, "Orbit size:", equiv, ". Number of constrains: ", &
nconstrain, ". Number of independent elemements:", nindep
                                    print *, ""
                                end if
                                do iaux = 1, 27
                                        do jaux = 1, 27
                                                kernel(ref3,iaux,jaux) = kern(iaux,jaux)
                                        end do
                                end do

                                n_indep_fc(ref3) = nindep
                                do iaux = 1, nindep
                                        indep_fc(ref3,iaux) = indep(iaux)-1 ! Giving it in python indexing
                                end do
                                if (nall3 == tot3) exit doii
                        end do dokk
                end do dojj
        end do doii
end subroutine

subroutine generate_rot4(rot_cart, Rot, nsym)

    double precision, dimension(nsym,3,3), intent(in) :: rot_cart
    double precision, dimension(24,nsym,81,81), intent(out) :: Rot
    
    integer, dimension(24,4) :: perms
    integer, dimension(4) :: cindex

    integer :: nsym
    integer :: iperm,isym,alphap,betap,gammap,thetap,indexp,alpha,beta,gamma,theta,index
    integer :: alphaperm, betaperm, gammaperm, thetaperm 

    perms = reshape([ &
        1,2,3,4,  1,2,4,3,  1,3,2,4,  1,3,4,2,  1,4,2,3,  1,4,3,2, & ! 1 at start
        2,1,3,4,  2,1,4,3,  2,3,1,4,  2,3,4,1,  2,4,1,3,  2,4,3,1, & ! 2 at start
        3,1,2,4,  3,1,4,2,  3,2,1,4,  3,2,4,1,  3,4,1,2,  3,4,2,1, & ! 3 at start
        4,1,2,3,  4,1,3,2,  4,2,1,3,  4,2,3,1,  4,3,1,2,  4,3,2,1 ], & ! 4 at start
        [24, 4], order=[2,1])

    do iperm = 1, 24
        do isym = 1, nsym
            do alphap = 1, 3
                do betap = 1, 3
                    do gammap = 1, 3
                        do thetap = 1, 3
                            indexp = 3*(3*(3*(alphap-1)+betap-1)+gammap-1)+thetap
                            do alpha = 1, 3
                                cindex(1) = alpha
                                do beta = 1, 3
                                    cindex(2) = beta
                                    do gamma = 1, 3
                                        cindex(3) = gamma
                                        do theta = 1, 3
                                            cindex(4) = theta
                                            index = 3*(3*(3*(alpha-1)+beta-1)+gamma-1)+theta
                                            alphaperm=cindex(perms(iperm,1))
                                            betaperm=cindex(perms(iperm,2))
                                            gammaperm=cindex(perms(iperm,3))
                                            thetaperm=cindex(perms(iperm,4))
                                            Rot(iperm,isym,indexp,index) = &
rot_cart(isym,alphap,alphaperm)*rot_cart(isym,betap,betaperm)*rot_cart(isym,gammap,gammaperm)*rot_cart(isym,thetap,thetaperm)
                                        end do
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
end subroutine generate_Rot4

subroutine recognize_quadruplet(nat, nat_sc, tot4, nref4, nsym, mappings, map_uc, nontrivial, M, verbose, &
        orbit4t, orbit4o, norbit, indep_fc, n_indep_fc, kernel, mapping_quadruplet)

        integer, intent(in) :: nat, nat_sc, nsym, tot4, nref4
        integer, dimension(:,:,:), intent(in) :: nontrivial
        double precision, dimension(:,:,:,:), intent(in) :: M
        integer, dimension(:,:), intent(in) :: mappings, map_uc
        logical, intent(in) :: verbose

        double precision, dimension(nref4,81,81), intent(out) :: kernel
        integer, dimension(nref4,81), intent(out) :: indep_fc 
        integer, dimension(nref4), intent(out) :: norbit, n_indep_fc
        integer, dimension(nref4,24*nsym,4), intent(out) :: orbit4t
        integer, dimension(nref4,24*nsym,2), intent(out) :: orbit4o
        integer, dimension(nat,nat_sc,nat_sc,nat_sc,3), intent(out) :: mapping_quadruplet

        double precision, allocatable, dimension(:,:) :: constrain_reduced

        integer, dimension(tot4,4) :: all4
        integer, dimension(24*nsym, 4) :: equilist
        integer, dimension(81) :: indep
        double precision, dimension(81,81) :: kern
        double precision, dimension(24*nsym*81,81) :: constrain

        integer, dimension(24,4) :: permutations
        integer, dimension(4) :: qplet, qplet_perm, qplet_sym
        integer :: ii, jj, kk, ll, nall4, equiv, nconstrain, iperm, isym, iaux, &
jaux, indexprime, mm, nindep, iw, ref4
        logical :: its_in_list

        character(len=100) :: filename


        kernel = 0
        orbit4t = 0
        orbit4o = 0
        mapping_quadruplet = 0

permutations = reshape([ &
        1,2,3,4,  1,2,4,3,  1,3,2,4,  1,3,4,2,  1,4,2,3,  1,4,3,2, & ! 1 at start
        2,1,3,4,  2,1,4,3,  2,3,1,4,  2,3,4,1,  2,4,1,3,  2,4,3,1, & ! 2 at start
        3,1,2,4,  3,1,4,2,  3,2,1,4,  3,2,4,1,  3,4,1,2,  3,4,2,1, & ! 3 at start
        4,1,2,3,  4,1,3,2,  4,2,1,3,  4,2,3,1,  4,3,1,2,  4,3,2,1 ], & ! 4 at start
        [24, 4], order=[2,1])

        ref4 = 0
        nall4 = 0
        doii : do ii = 0, nat-1
                dojj : do jj = 0, nat_sc-1
                        dokk: do kk = 0, nat_sc-1
                                doll: do ll = 0, nat_sc-1
                                        qplet = [ii, jj, kk, ll]
                                        call qplet_in_list(qplet, all4, nall4, its_in_list)
                                        if (its_in_list) cycle doll
                                        ref4 = ref4 + 1
                                        equiv = 0
                                        nconstrain = 0
                                        constrain = 0
                                        do iperm = 1, 24
                                                qplet_perm(1) = qplet(permutations(iperm,1))
                                                qplet_perm(2) = qplet(permutations(iperm,2))
                                                qplet_perm(3) = qplet(permutations(iperm,3))
                                                qplet_perm(4) = qplet(permutations(iperm,4))
                                                do isym = 1, nsym
                                                        qplet_sym = &
[mappings(qplet_perm(1)+1,isym), mappings(qplet_perm(2)+1,isym), mappings(qplet_perm(3)+1,isym), mappings(qplet_perm(4)+1,isym)]
                                                        if (qplet_sym(1) >= nat) then
                                                                qplet_sym(2) = map_uc(qplet_sym(1)+1, qplet_sym(2)+1)
                                                                qplet_sym(3) = map_uc(qplet_sym(1)+1, qplet_sym(3)+1)
                                                                qplet_sym(4) = map_uc(qplet_sym(1)+1, qplet_sym(4)+1)
                                                                qplet_sym(1) = map_uc(qplet_sym(1)+1, qplet_sym(1)+1)
                                                        end if
                                                        call qplet_in_list(qplet_sym, equilist, equiv, its_in_list)
                                                        if (((iperm==1) .and. (isym==1)) .or. (.not. its_in_list)) then
                                                                equilist(equiv+1,:) = qplet_sym
                                                                orbit4t(ref4,equiv+1,:) = qplet_sym
                                                                all4(nall4+1,:) = qplet_sym
                                                                orbit4o(ref4,equiv+1,:) = [iperm-1, isym-1]
                                                                mapping_quadruplet( &
qplet_sym(1)+1,qplet_sym(2)+1,qplet_sym(3)+1,qplet_sym(4)+1,:) = [ref4-1,iperm-1,isym-1]
                                                                equiv = equiv + 1
                                                                nall4 = nall4 + 1
                                                        end if
                                                        if (all(qplet == qplet_sym)) then
                                                                do indexprime = 1, 81
                                                                        if (nontrivial(iperm,isym,indexprime)==1) then
                                                                                do mm = 1, 81
                                                                                        constrain(nconstrain+1, mm) = & 
                                                                                M(iperm,isym,indexprime,mm)
                                                                                end do
                                                                                nconstrain = nconstrain + 1
                                                                        end if
                                                                end do
                                                        end if
                                                end do
                                        end do
                                        norbit(ref4) = equiv
                                        allocate(constrain_reduced(max(nconstrain,81),81))
                                        constrain_reduced(:,:) = 0
                                        do iaux = 1, nconstrain
                                                do jaux =1, 81
                                                        constrain_reduced(iaux,jaux) = constrain(iaux,jaux)
                                                end do
                                        end do
                                        call gauss_jordan(&
        constrain_reduced, max(nconstrain,81), 81, nconstrain, kern, indep, nindep)
                                        deallocate(constrain_reduced)
                                        if (verbose) then
                                            print*, &
        "Reference quadruplet: (", ii, ",", jj, ",", kk, ",", ll, ")."
                                            print*, "Orbit size:", equiv, ". Number of constrains: ", nconstrain, &
        ". Number of independet elements:", nindep
                                            print*, ""
                                        end if
                                        do iaux = 1, 81
                                                do jaux = 1, 81
                                                        kernel(ref4,iaux,jaux) = kern(iaux,jaux)
                                                end do
                                        end do

                                        n_indep_fc(ref4) = nindep
                                        do iaux = 1, nindep
                                                indep_fc(ref4,iaux) = indep(iaux)-1 ! Giving it in python indexing
                                        end do
                                        if (nall4 == tot4) exit doii
                                end do doll
                        end do dokk
                end do dojj
        end do doii
end subroutine   

subroutine triplet_in_list(triplet, all3, nall3, its_in_list)
        integer, dimension(3), intent(in) :: triplet
        integer, dimension(:,:), intent(in) :: all3
        integer, intent(in) :: nall3

        logical, intent(out) :: its_in_list
        
        integer :: i
        its_in_list=.false.
        do i = 1, nall3
                if (all(triplet == all3(i,:))) then
                        its_in_list=.true.
                        return
                end if
        end do
end subroutine

subroutine qplet_in_list(qplet, all4, nall4, its_in_list)
        integer, dimension(4), intent(in) :: qplet
        integer, dimension(:,:), intent(in) :: all4
        integer, intent(in) :: nall4

        logical, intent(out) :: its_in_list
        
        integer :: i
        its_in_list=.false.
        do i = 1, nall4
                if (all(qplet == all4(i,:))) then
                        its_in_list=.true.
                        return
                end if
        end do
end subroutine

subroutine is_qcons(qsum, its_cons)
    double precision, dimension(3), intent(in) :: qsum
    logical, intent(out) :: its_cons

    double precision :: q1, q2, q3
    double precision, parameter :: eps=1e-3

    q1 = modulo(abs(qsum(1))+eps,1.0d0)
    q2 = modulo(abs(qsum(2))+eps,1.0d0)
    q3 = modulo(abs(qsum(3))+eps,1.0d0)
    if ((q1<2*eps) .and. (q2<2*eps) .and. (q3<2*eps)) then
        its_cons = .true.
    else
        its_cons = .false.
    end if
end subroutine

subroutine gauss_jordan(a, lda, n, nconstrain, b, independent, n_indep)
    implicit none

    integer, intent(in) :: lda, n, nconstrain
    real(8), intent(inout) :: a(lda, n)
    real(8), intent(out) :: b(n, n)
    integer, intent(out) :: independent(n)
    integer, intent(out) :: n_indep

    real(8), parameter :: EPS = 1d-10
    integer :: dependent(n)
    integer :: i, j, k, i_row, pivot_row, n_dep
    real(8) :: tmp, pivot_val

    ! Initialize
    b = 0.0d0
    n_dep = 0
    n_indep = 0
    i_row = 1 

    do k = 1, n
        pivot_row = 0
        
        ! Only look for a pivot if we have rows left within the 13
        if (i_row <= nconstrain) then
            pivot_row = i_row
            ! Find max in column k (rows i_row to 13)
            do i = i_row + 1, nconstrain
                if (abs(a(i, k)) > abs(a(pivot_row, k))) pivot_row = i
            end do

            if (abs(a(pivot_row, k)) > EPS) then
                ! SUCCESS: Found a pivot
                n_dep = n_dep + 1
                dependent(n_dep) = k

                ! 1. Swap rows to move pivot to i_row
                if (pivot_row /= i_row) then
                    do j = 1, n
                        tmp = a(i_row, j)
                        a(i_row, j) = a(pivot_row, j)
                        a(pivot_row, j) = tmp
                    end do
                end if

                ! 2. Normalize the pivot row
                pivot_val = a(i_row, k)
                do j = 1, n
                    a(i_row, j) = a(i_row, j) / pivot_val
                end do
                ! Force exact 1.0 to avoid precision drift
                a(i_row, k) = 1.0d0

                ! 3. Eliminate other rows (only within the 13 logical rows)
                do i = 1, nconstrain
                    if (i /= i_row) then
                        tmp = a(i, k)
                        if (abs(tmp) > EPS) then
                            do j = 1, n
                                a(i, j) = a(i, j) - tmp * a(i_row, j)
                            end do
                            a(i, k) = 0.0d0 ! Force exact 0.0
                        end if
                    end if
                end do
                
                i_row = i_row + 1
                cycle 
            end if
        end if

        ! No pivot found in this column -> Independent
        n_indep = n_indep + 1
        independent(n_indep) = k
    end do

    ! 4. Build Null Space Basis B
    ! Only map if we actually found pivots
    !if (n_dep > 0) then
    do j = 1, n_indep
            do i = 1, n_dep
            ! Row 'i' of RREF matrix 'a' corresponds to the 'i-th' dependent column
            b(dependent(i), j) = -a(i, independent(j))
            end do
            b(independent(j), j) = 1.0d0
    end do
    !end if

end subroutine gauss_jordan

end module module_hess

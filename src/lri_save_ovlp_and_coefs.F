!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2020 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Save LRI overlap matrix  and expansion coefficients
!>        lri : local resolution of the identity
!> \par History
!>      created Andrea Grisafi [04.2021]
!> \authors Andrea Grisafi 
! **************************************************************************************************
MODULE lri_save_ovlp_and_coefs
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind,&
                                              get_atomic_kind_set
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type,&
                                              init_orb_basis_set
   USE cell_types,                      ONLY: cell_type
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_generate_filename,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_get_block_p,&
                                              dbcsr_p_type,&
                                              dbcsr_type
   USE generic_os_integrals,            ONLY: int_overlap_aabb_os
   USE generic_shg_integrals,           ONLY: int_overlap_ab_shg_low, &
                                              lri_precalc_angular_shg_part
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE lri_environment_init,            ONLY: lri_basis_init
   USE lri_environment_methods,         ONLY: calculate_avec_lri,&
                                              calculate_lri_integrals
   USE lri_environment_types,           ONLY: allocate_lri_ints_rho,&
                                              deallocate_lri_ints_rho,&
                                              lri_density_type,&
                                              lri_environment_type,&
                                              lri_int_rho_type,&
                                              lri_int_type,&
                                              lri_list_type,&
                                              lri_rhoab_type,&
                                              lri_kind_type
   USE lri_save_ovlp_and_coefs_types,   ONLY: create_lri_save,&
                                              deallocate_lri_save,&
                                              lri_save_type

   USE memory_utilities,                ONLY: reallocate
   USE message_passing,                 ONLY: mp_sum
   USE particle_types,                  ONLY: particle_type
   USE powell,                          ONLY: opt_state_type,&
                                              powell_optimize
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type,&
                                              set_qs_env
   USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                              neighbor_list_iterate,&
                                              neighbor_list_iterator_create,&
                                              neighbor_list_iterator_p_type,&
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! **************************************************************************************************

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'lri_save_ovlp_and_coefs'

   PUBLIC :: save_ovlp_and_coefs

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief save lri ovlp and coefs 
!> \param qs_env qs environment
! **************************************************************************************************
   SUBROUTINE save_ovlp_and_coefs(qs_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env

      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: pmatrix
      TYPE(lri_density_type), POINTER                    :: lri_density
      TYPE(lri_environment_type), POINTER                :: lri_env
      TYPE(opt_state_type)                               :: opt_state
      TYPE(qs_rho_type), POINTER                         :: rho_struct
      TYPE(section_vals_type), POINTER                   :: dft_section, input
      TYPE(lri_save_type), POINTER                       :: lri_save
      TYPE(section_vals_type), POINTER                   :: lri_save_section
      TYPE(lri_kind_type), DIMENSION(:), POINTER         :: lri_coef
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind
      TYPE(gto_basis_set_type), POINTER                  :: fbas,fbasa,fbasb,obasa,obasb
      REAL(KIND=dp), DIMENSION(500, 500)                 :: overlap 
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: ovlp,selfovlp 
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: scon_fba, scon_fbb
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: dummy1, Waux_mat
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :)  :: dummy2, dWaux_mat
      REAL(KIND=dp), DIMENSION(3)                        :: rab
      REAL(KIND=dp)                                      :: dab 
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: soo_list

      !INTEGER, ALLOCATABLE, DIMENSION(:,:)               :: irep
      INTEGER                                            :: atom_a, atom_b, iatom, iunit, &
                                                            ikind, ilist, ispin, jatom, jkind, &
                                                            jneighbor, nat, natom, nfa, nfb, &
                                                            nkind, nspin, nsgf, icoeff, ifa, ifb, &
                                                            mepos, nneighbor, nlist, nthread, &
                                                            iatom_old, jatom_old

      CHARACTER*16                                       :: ikind_str, iatom_str, jatom_str, rep_str
      CHARACTER*64                                       :: filename 

      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, input=input, &
                      lri_env=lri_env, lri_density=lri_density, nkind=nkind, &
                      para_env=para_env, rho=rho_struct)

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, natom=nat)
      
      iatom_old = 1
      jatom_old = 1
      overlap = 0.d0
      mepos = 0
      nthread = 1

      NULLIFY(nl_iterator,soo_list,fbasa,fbasb,scon_fba,scon_fbb,obasa,obasb)

      IF (ASSOCIATED(lri_env%soo_list)) THEN

         soo_list => lri_env%soo_list
         CALL neighbor_list_iterator_create(nl_iterator, soo_list, nthread=nthread)

         DO WHILE (neighbor_list_iterate(nl_iterator, mepos) == 0)
         
            CALL get_iterator_info(nl_iterator, mepos=mepos, ikind=ikind, jkind=jkind, &
                                   nlist=nlist, ilist=ilist, nnode=nneighbor, inode=jneighbor, &
                                   iatom=iatom, jatom=jatom, r=rab)
           
            IF (iatom.ne.iatom_old .or. jatom.ne.jatom_old) THEN
               overlap = 0.d0
            ENDIF 
 
            obasa => lri_env%orb_basis(ikind)%gto_basis_set
            obasb => lri_env%orb_basis(jkind)%gto_basis_set
            fbasa => lri_env%ri_basis(ikind)%gto_basis_set
            fbasb => lri_env%ri_basis(jkind)%gto_basis_set
            scon_fba => lri_env%bas_prop(ikind)%scon_ri
            scon_fbb => lri_env%bas_prop(jkind)%scon_ri

            CALL lri_precalc_angular_shg_part(obasa, obasb, fbasa, fbasb, &
                  rab, Waux_mat, dWaux_mat,.FALSE.)
           
            nfa = fbasa%nsgf
            nfb = fbasb%nsgf
            ALLOCATE(ovlp(nfa,nfb))
            IF (ALLOCATED(ovlp)) THEN
               CALL int_overlap_ab_shg_low(ovlp, dummy1, rab, fbasa, fbasb, scon_fba, scon_fbb, &
                                        Waux_mat, dWaux_mat, .TRUE., .FALSE., contraction_high=.FALSE.)
            ENDIF           
            overlap(:nfa,:nfb) = overlap(:nfa,:nfb) + ovlp(:nfa,:nfb)
            DEALLOCATE(ovlp)

            write(iatom_str,*) iatom
            write(jatom_str,*) jatom
            filename = "overlaps/overlap_"//trim(adjustl(iatom_str))//"-"&
                                          //trim(adjustl(jatom_str))//".dat"
                                          !//trim(adjustl(jatom_str))//"_rep"&
                                          !//trim(adjustl(rep_str))//".dat"
            CALL open_file(filename, unit_number=iunit, file_status="UNKNOWN",file_action="WRITE")
            do ifa=1,nfa
              do ifb=1,nfb
                write(iunit,*) overlap(ifa, ifb)
              enddo
            enddo
            CALL close_file(iunit)

            iatom_old = iatom
            jatom_old = jatom
 
          ENDDO

      ENDIF

      !DEALLOCATE(irep)     
 
      ! SAVE LRI COEFFICIENTS
      ! -----------------------------------------------------------
      ! get lri coefs: closed-shell system assumed (ispin=1)
      lri_coef => lri_density%lri_coefs(1)%lri_kinds
      ! for each atomic type
      DO ikind = 1, nkind
         write(ikind_str,*) ikind
         ! get the atomic kind
         atomic_kind => atomic_kind_set(ikind)
         ! RI auxiliary basis
         fbas => lri_env%ri_basis(ikind)%gto_basis_set
         ! number of RI spherical Gaussian functions
         nsgf = fbas%nsgf
         ! get the number of atoms associated with atomic kind
         CALL get_atomic_kind(atomic_kind=atomic_kind, natom=natom)
         ! for each atom of that type
         DO iatom = 1, natom
            write(iatom_str,*) iatom
            ! open file to save coefficients
            filename="coefficients/coefs_type"//trim(adjustl(ikind_str))//"_atom"&
                                              //trim(adjustl(iatom_str))//".dat"
            CALL open_file(filename, unit_number=iunit, file_status="UNKNOWN",file_action="WRITE")
            ! for each coef
            DO icoeff=1,nsgf
               ! save coef
               write(iunit,*) lri_coef(ikind)%acoef(iatom, icoeff)
            END DO !for each coef
            ! close file
            CALL close_file(iunit)
         END DO ! for each atom of that type
      END DO ! for each atomic type


   END SUBROUTINE save_ovlp_and_coefs 

END MODULE lri_save_ovlp_and_coefs

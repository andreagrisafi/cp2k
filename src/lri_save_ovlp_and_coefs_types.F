!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2020 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief sets the environment for saving the lri overlap and coefficients 
!>        lri : local resolution of the identity
!> \par History
!>      created Andrea Grisafi [04.2021]
!> \authors Andrea Grisafi
! **************************************************************************************************
MODULE lri_save_ovlp_and_coefs_types

   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: pi
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'lri_save_ovlp_and_coefs'
   PUBLIC :: lri_save_type
   PUBLIC :: create_lri_save, deallocate_lri_save

! **************************************************************************************************
   TYPE lri_save_type
      LOGICAL                                            :: use_save_lri
   END TYPE lri_save_type

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief creates lri_save
!> \param lri_save environment
! **************************************************************************************************
   SUBROUTINE create_lri_save(lri_save)

      TYPE(lri_save_type), POINTER                        :: lri_save

      ALLOCATE (lri_save)

      lri_save%use_save_lri = .TRUE.

   END SUBROUTINE create_lri_save

! **************************************************************************************************
!> \brief deallocates lri_save
!> \param lri_save environment
! **************************************************************************************************
   SUBROUTINE deallocate_lri_save(lri_save)

      TYPE(lri_save_type), POINTER                        :: lri_save

      INTEGER                                            :: i

      IF (ASSOCIATED(lri_save)) THEN
         DEALLOCATE (lri_save)
      ENDIF
   END SUBROUTINE deallocate_lri_save

END MODULE lri_save_ovlp_and_coefs_types

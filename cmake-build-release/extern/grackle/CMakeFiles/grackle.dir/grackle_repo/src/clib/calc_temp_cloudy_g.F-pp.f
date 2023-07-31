# 1 "/Users/noah/Documents/SPH/4 Jan 22/SPH-EXA-fork/extern/grackle/grackle_repo/src/clib/calc_temp_cloudy_g.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/Users/noah/Documents/SPH/4 Jan 22/SPH-EXA-fork/extern/grackle/grackle_repo/src/clib/calc_temp_cloudy_g.F"

# 1 "/Users/noah/Documents/SPH/4 Jan 22/SPH-EXA-fork/extern/grackle/grackle_repo/src/clib/phys_const.def" 1

# 19 "/Users/noah/Documents/SPH/4 Jan 22/SPH-EXA-fork/extern/grackle/grackle_repo/src/clib/phys_const.def"

# 33 "/Users/noah/Documents/SPH/4 Jan 22/SPH-EXA-fork/extern/grackle/grackle_repo/src/clib/phys_const.def"



# 2 "/Users/noah/Documents/SPH/4 Jan 22/SPH-EXA-fork/extern/grackle/grackle_repo/src/clib/calc_temp_cloudy_g.F" 2

!=======================================================================
!//////////////////  SUBROUTINE CALC_TEMP_CLOUDY_G  \\\\\\\\\\\\\\\\\\\

      subroutine calc_temp_cloudy_g(d, e, metal, temperature,
     &                in, jn, kn, iexpand, imetal,
     &                is, js, ks, ie, je, ke,
     &                aye, temstart, temend, 
     &                utem, uxyz, uaye, urho, utim,
     &                gamma, fh, 
     &                priGridRank, priGridDim,
     &                priPar1, priPar2, priPar3,
     &                priDataSize, priMMW)

!
!  CALCULATE TEMPERATURE USING MMW FROM CLOUDY TABLE
!
!  written by: Britton Smith
!  date:       May, 2015
!
!  PURPOSE:
!    Calculate the temperature field using tabulated mu
!
!  INPUTS:
!    in,jn,kn - dimensions of 3D fields
!
!    d        - total density field
!    e        - internal energy field
!    metal    - metal density field
!    temperature - temperature field
!
!    is,ie    - start and end indices of active region (zero based)
!    iexpand  - comoving coordinates flag (0 = off, 1 = on)
!    imetal   - flag if metal field is active (0 = no, 1 = yes)
!
!    fh       - Hydrogen mass fraction (typically 0.76)
!    aye      - expansion factor (in code units)
!
!    utim     - time units (i.e. code units to CGS conversion factor)
!    uaye     - expansion factor conversion factor (uaye = 1/(1+zinit))
!    urho     - density units
!    uxyz     - length units
!    utem     - temperature(-like) units
!
!    temstart, temend - start and end of temperature range for rate table
!
!    priGridRank - rank of cloudy primordial cooling data grid
!    priGridDim  - array containing dimensions of cloudy primordial data
!    priPar1, priPar2, priPar3 - arrays containing primordial grid parameter values
!    priDataSize - total size of flattened 1D primordial cooling data array
!    priMMW      - primordial mmw data
!
!  OUTPUTS:
!    update temperature array
!
!  PARAMETERS:
!    mh      - H mass in cgs units
!
!-----------------------------------------------------------------------

      implicit NONE

# 1 "/Users/noah/Documents/SPH/4 Jan 22/SPH-EXA-fork/extern/grackle/grackle_repo/src/clib/grackle_fortran_types.def" 1
!=======================================================================
!
!
! Grackle fortran variable types
!
!
! Copyright (c) 2013, Enzo/Grackle Development Team.
!
! Distributed under the terms of the Enzo Public Licence.
!
! The full license is in the file LICENSE, distributed with this 
! software.
!=======================================================================












      integer, parameter :: RKIND=8





      integer, parameter :: DKIND=8
      integer, parameter :: DIKIND=8
# 64 "/Users/noah/Documents/SPH/4 Jan 22/SPH-EXA-fork/extern/grackle/grackle_repo/src/clib/calc_temp_cloudy_g.F" 2




!  General Arguments

      integer in, jn, kn, is, js, ks, ie, je, ke, 
     &        iexpand, imetal
      real*8  aye, temstart, temend, gamma,
     &        utim, uxyz, uaye, urho, utem, fh

!  Density, energy and velocity fields fields

      real*8  d(in,jn,kn), e(in,jn,kn), metal(in,jn,kn),
     &     temperature(in,jn,kn)

!  Cloudy cooling data

      integer*8 priGridRank, priDataSize,
     &     priGridDim(priGridRank)
      real*8 priPar1(priGridDim(1)), priPar2(priGridDim(2)), 
     &     priPar3(priGridDim(3)), priMMW(priDataSize)

!  Parameters

      real*8 mh
      parameter (mh = 1.67262171d-24)

!  Locals

      integer i, j, k
      integer t, dj, dk
      real*8 dom, zr
      real*8 dbase1, tbase1, xbase1

!  row temporaries

      real*8 tgas(in), rhoH(in), mmw(in)

!  Iteration mask

      logical itmask(in)
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
!=======================================================================
      
!     Set units

      dom      = urho*(aye**3)/mh
      tbase1   = utim
      xbase1   = uxyz/(aye*uaye)    ! uxyz is [x]*a      = [x]*[a]*a'        '
      dbase1   = urho*(aye*uaye)**3 ! urho is [dens]/a^3 = [dens]/([a]*a')^3 '
      zr       = 1._DKIND/(aye*uaye) - 1._DKIND

!  Convert densities from comoving to proper

      if (iexpand .eq. 1) then

         call scale_fields_table_g(d, metal,
     &                  is, ie, js, je, ks, ke,
     &                  in, jn, kn, imetal, aye**(-3))

      endif

!  Loop over zones, and do an entire i-column in one go

      dk = ke - ks + 1
      dj = je - js + 1

! parallelize the k and j loops with OpenMP
! flat j and k loops for better parallelism






      do t = 0, dk*dj-1
        k = t/dj      + ks+1
        j = mod(t,dj) + js+1

!     Initialize iteration mask to true for all cells.

         do i = is+1, ie+1
            itmask(i) = .true.

            rhoH(i) = fh * d(i,j,k)
         enddo

!     Calculate temperature and mean molecular weight

         call calc_temp1d_cloudy_g(d, metal, e, rhoH,
     &        in, jn, kn, is, ie, j, k,
     &        tgas, mmw, dom, zr, 
     &        temstart, temend,
     &        gamma, utem, imetal,
     &        priGridRank, priGridDim,
     &        priPar1, priPar2, priPar3,
     &        priDataSize, priMMW,
     &        itmask)

!     Copy slice values into field array

         do i = is+1, ie+1
            temperature(i,j,k) = tgas(i)
         enddo

      enddo




!     Convert densities back to comoving from proper

      if (iexpand .eq. 1) then

         call scale_fields_table_g(d, metal,
     &                  is, ie, js, je, ks, ke,
     &                  in, jn, kn, imetal, aye**3)

      endif

      return
      end

c -----------------------------------------------------------
!   This routine scales the density fields from comoving to
!     proper densities (and back again).

      subroutine scale_fields_table_g(d, metal,
     &                        is, ie, js, je, ks, ke,
     &                        in, jn, kn, imetal, factor)
c -------------------------------------------------------------------

      implicit NONE

# 1 "/Users/noah/Documents/SPH/4 Jan 22/SPH-EXA-fork/extern/grackle/grackle_repo/src/clib/grackle_fortran_types.def" 1
!=======================================================================
!
!
! Grackle fortran variable types
!
!
! Copyright (c) 2013, Enzo/Grackle Development Team.
!
! Distributed under the terms of the Enzo Public Licence.
!
! The full license is in the file LICENSE, distributed with this 
! software.
!=======================================================================












      integer, parameter :: RKIND=8





      integer, parameter :: DKIND=8
      integer, parameter :: DIKIND=8
# 200 "/Users/noah/Documents/SPH/4 Jan 22/SPH-EXA-fork/extern/grackle/grackle_repo/src/clib/calc_temp_cloudy_g.F" 2

!     Arguments

      integer in, jn, kn, is, ie, js, je, ks, ke, imetal
      real*8  d(in,jn,kn), metal(in,jn,kn)
      real*8 factor

!     locals

      integer i, j, k

!     Multiply all fields by factor (1/a^3 or a^3)

      do k = ks+1, ke+1
         do j = js+1, je+1
            do i = is+1, ie+1
               d(i,j,k)        = d(i,j,k) * factor
            enddo
         enddo
      enddo

      if (imetal .eq. 1) then
         do k = ks+1, ke+1
            do j = js+1, je+1
               do i = is+1, ie+1
                  metal(i,j,k) = metal(i,j,k) * factor
               enddo
            enddo
         enddo
      endif

      return
      end

MODULE FLUE_philox_helpers
   USE FLUE_constants, ONLY : WP
   !!
   !! Helper utilities for the Philox-based heatbath implementation.
   !!
   !! This module provides:
   !!   * a uniform 64-bit integer kind alias (C64),
   !!   * deterministic derivation of a stage-specific key from a master key,
   !!   * a linear site index for mapping a lattice site to a Philox counter,
   !!   * colouring rules for Wilson (2-colour) and Symanzik (4-colour) updates.
  !!
  USE Philox, ONLY: C64
   IMPLICIT NONE(TYPE, EXTERNAL)
   PRIVATE

   PUBLIC :: derive_stage_key
   PUBLIC :: site_linear_index
   PUBLIC :: site_colour
   PUBLIC :: ncolours_for_action

CONTAINS

   SUBROUTINE derive_stage_key(master_key, sweep_id, stage_tag, mu, colour, key)
     !! Inputs:
     !!   master_key : 128-bit base key (2 x 64-bit words)
     !!   sweep_id   : trajectory / sweep index
     !!   stage_tag  : algorithm stage identifier
     !!   mu         : link direction
     !!   colour      : checkerboard / colour subset
     !!
     !! Output:
     !!   key        : derived 128-bit Philox key

     INTEGER(C64), INTENT(IN)  :: master_key(2)
     INTEGER,      INTENT(IN)  :: sweep_id, stage_tag, mu, colour
     INTEGER(C64), INTENT(OUT) :: key(2)
     INTEGER(C64) :: t1, t2
     !-----------------------------------------------------------
     ! Construct t1 by packing (sweep_id, mu, colour) into a 64-bit integer
     !
     ! shiftl(x, n) ≡ x * 2^n
     !
     ! So:
     !   shiftl(int(sweep_id), 32) ≡ sweep_id * 2^32
     !   shiftl(int(mu),       16) ≡ mu       * 2^16
     !
     ! Therefore:
     !   t1 = sweep_id * 2^32 + mu * 2^16 + colour
     !
     ! Bit layout:
     !   [ sweep_id (upper 32 bits) | mu (next 16 bits) | colour (low 16 bits) ]
     !---------------------------------------------------------
     t1 = shiftl(int(sweep_id, C64), 32) &   ! = sweep_id * 2^32
          + shiftl(int(mu, C64), 16) &   ! = mu       * 2^16
          + int(colour, C64)                   ! = colour
     !-------------------------------------------------------------
     ! Construct t2 by packing (stage_tag, mu, colour) into a different layout
     !
     !   shiftl(int(stage_tag), 40) ≡ stage_tag * 2^40
     !   shiftl(int(mu),        20) ≡ mu        * 2^20
     !
     ! Therefore:
     !   t2 = stage_tag * 2^40 + mu * 2^20 + colour
     !
     ! Bit layout:
     !   [ stage_tag (upper bits) | mu (next 20 bits) | colour (low 20 bits) ]
     !
     ! Note: different shifts from t1 intentionally decorrelate fields.
     !----------------------------------------------------------
     t2 = shiftl(int(stage_tag, C64), 40) &  ! = stage_tag * 2^40
          + shiftl(int(mu, C64), 20) &  ! = mu        * 2^20
          + int(colour, C64)                   ! = colour
     !---------------------------------------------------------
     ! Mix into final 128-bit key using XOR
     !
     ! key(1):
     !   key(1) = master_key(1) XOR t1
     !
     ! key(2):
     !   shiftl(t1, 7) ≡ t1 * 2^7 = t1 * 128
     !
     !   key(2) = master_key(2) XOR [ t2 XOR (t1 * 128) ]
     !
     ! The shift by 7 spreads t1 bits into different positions before mixing,
     ! reducing linear correlation between key(1) and key(2).
     !---------------------------------------------------------
     key(1) = ieor(master_key(1), t1)
     key(2) = ieor(master_key(2), &
          ieor(t2, shiftl(t1, 7)))   ! shiftl(t1,7) = t1 * 2^7 = t1 * 128
   END SUBROUTINE derive_stage_key


   PURE INTEGER FUNCTION ncolours_for_action(use_symanzik) RESULT(ncolours)
      !!
      !! Return the number of independent colour classes needed for the action.
      !!
      !! Wilson   -> 2 colours
      !! Symanzik -> 4 colours
      !!
      !! The Symanzik rectangles couple same-checkerboard sites, so 2 colours are
      !! not sufficient there.
      !!
      LOGICAL, INTENT(IN) :: use_symanzik

      IF (use_symanzik) THEN
         ncolours = 4
      ELSE
         ncolours = 2
      END IF
   END FUNCTION ncolours_for_action


   PURE INTEGER FUNCTION site_colour(coord, mu, use_symanzik) RESULT(colour)
      !!
      !! Compute the colour class of a site for a fixed link direction mu.
      !!
      !! coord = [it, ix, iy, iz]
      !!
      !! Wilson:
      !!   Use the ordinary checkerboard parity:
      !!      colour = mod(it + ix + iy + iz, 2)
      !!
      !! Symanzik:
      !!   Use a 4-colour partition that separates the displacements which appear
      !!   in the rectangular loops:
      !!      colour = 2*mod(coord(mu),2) + mod(sum(other coords),2)
      !!
      INTEGER, INTENT(IN) :: coord(4), mu
      LOGICAL, INTENT(IN) :: use_symanzik
      INTEGER :: transverse_parity

      IF (use_symanzik) THEN
         transverse_parity = mod(sum(coord) - coord(mu), 2)
         colour = 2 * mod(coord(mu), 2) + transverse_parity
      ELSE
         colour = mod(sum(coord), 2)
      END IF
   END FUNCTION site_colour

   PURE INTEGER(C64) FUNCTION site_linear_index(coord, dims) RESULT(idx)
     !!
     !! Convert a 4D lattice coordinate to a unique 1-based linear site index.
     !!
     !! coord = [it, ix, iy, iz]
     !! dims  = [nt, nx, ny, nz]
     !!
     !! This linear index is used as counter(1) for the Philox stream so each
     !! site has its own independent random-number namespace.
     !!
     INTEGER, INTENT(IN) :: coord(4), dims(4)
     idx = int((((coord(1) - 1) * dims(2) + (coord(2) - 1)) * dims(3) + &
          (coord(3) - 1)) * dims(4) + (coord(4) - 1), C64) + 1_C64
   END FUNCTION site_linear_index

END MODULE FLUE_philox_helpers

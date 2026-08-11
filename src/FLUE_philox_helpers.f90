module FLUE_philox_helpers
  !<
  !< Helper utilities for the Philox-based heatbath implementation.
  !<
  !< This module provides:
  !<   * a uniform 64-bit integer kind alias (C64),
  !<   * deterministic derivation of a stage-specific key from a master key,
  !<   * a linear site index for mapping a lattice site to a Philox counter,
  !<   * colouring rules for Wilson (2-colour) and Symanzik (4-colour) updates.
  !<
  use FLUE_constants, only: WP
   use Philox, only: C64
   implicit none(type, external)
   private

   public :: derive_stage_key
   public :: site_linear_index
   public :: site_colour
   public :: ncolours_for_action

contains

  subroutine derive_stage_key(master_key, sweep_id, stage_tag, mu, colour, key)
    !< Create a separated key for each stage
    !<
     !< Construct t1 by packing (sweep_id, mu, colour) into a 64-bit integer
     !<
     !< shiftl(x, n) ≡ x * 2^n
     !<
     !< So:
     !<   shiftl(int(sweep_id), 32) ≡ sweep_id * 2^32
     !<   shiftl(int(mu),       16) ≡ mu       * 2^16
     !<
     !< Therefore:
     !<   t1 = sweep_id * 2^32 + mu * 2^16 + colour
     !<
     !< Bit layout:
     !<   [ sweep_id (upper 32 bits) | mu (next 16 bits) | colour (low 16 bits) ]
    !<
     !< Construct t2 by packing (stage_tag, mu, colour) into a different layout
     !<
     !<   shiftl(int(stage_tag), 40) ≡ stage_tag * 2^40
     !<   shiftl(int(mu),        20) ≡ mu        * 2^20
     !<
     !< Therefore:
     !<   t2 = stage_tag * 2^40 + mu * 2^20 + colour
     !<
     !< Bit layout:
     !<   [ stage_tag (upper bits) | mu (next 20 bits) | colour (low 20 bits) ]
     !<
     !< Note: different shifts from t1 intentionally decorrelate fields.
    !<
     !< Mix into final 128-bit key using XOR
     !<
     !< key(1):
     !<   key(1) = master_key(1) XOR t1
     !<
     !< key(2):
     !<   shiftl(t1, 7) ≡ t1 * 2^7 = t1 * 128
     !<
     !<   key(2) = master_key(2) XOR [ t2 XOR (t1 * 128) ]
     !<
     !< The shift by 7 spreads t1 bits into different positions before mixing,
     !< reducing linear correlation between key(1) and key(2).
     integer(kind=C64), intent(IN), dimension(2) :: master_key
     !< 128-bit base key (2 x 64-bit words)
     integer, intent(IN) :: sweep_id
     !< trajectory / sweep index
     integer, intent(IN) :: stage_tag
     !< algorithm stage identifier
     integer, intent(IN) :: mu
     !< link direction
     integer, intent(IN) :: colour
     !< checkerboard subset
     integer(kind=C64), intent(OUT), dimension(2) :: key
     !< derived 128-bit Philox key
     integer(kind=C64) :: t1, t2
     !< sweep_id/mu/colour & stage_tag/mu/colour keys to be mixed
      t1 = SHIFTL(INT(sweep_id, C64), 32) &   ! = sweep_id * 2^32
           + SHIFTL(INT(mu, C64), 16) &   ! = mu       * 2^16
           + INT(colour, C64)                   ! = colour

      t2 = SHIFTL(INT(stage_tag, C64), 40) &  ! = stage_tag * 2^40
           + SHIFTL(INT(mu, C64), 20) &  ! = mu        * 2^20
           + INT(colour, C64)                   ! = colour

      key(1) = IEOR(master_key(1), t1)
      key(2) = IEOR(master_key(2), &
                    IEOR(t2, SHIFTL(t1, 7)))   ! shiftl(t1,7) = t1 * 2^7 = t1 * 128
   end subroutine derive_stage_key

   pure function ncolours_for_action(use_symanzik) result(ncolours)
      !<
      !< Return the number of independent colour classes needed for the action.
      !<
      !< Wilson   -> 2 colours
      !< Symanzik -> 4 colours
      !<
      !< The Symanzik rectangles couple same-checkerboard sites, so 2 colours are
      !< not sufficient there.
      !<
     logical, intent(IN) :: use_symanzik
     !< Whether we use a symanzik improved action or not
     integer :: ncolours
     !< The size of the checkerboard required
      if (use_symanzik) then
         ncolours = 4
      else
         ncolours = 2
      end if
   end function ncolours_for_action

   pure function site_colour(coord, mu, use_symanzik) result(colour)
      !<
      !< Compute the colour class of a site for a fixed link direction mu.
      !<
      !< coord = [it, ix, iy, iz]
      !<
      !< Wilson:
      !<   Use the ordinary checkerboard parity:
      !<      colour = mod(it + ix + iy + iz, 2)
      !<
      !< Symanzik:
      !<   Use a 4-colour partition that separates the displacements which appear
      !<   in the rectangular loops:
      !<      colour = 2*mod(coord(mu),2) + mod(sum(other coords),2)
      !<
     integer, intent(IN), dimension(4) :: coord
     !< coordinate of the link
     integer, intent(IN) :: mu
     !< integer direction of link
     logical, intent(IN) :: use_symanzik
     !< The actioin type
     integer :: colour
     !< checkerboard index
     integer :: transverse_parity
     !< for calculating checkerboard
     if (use_symanzik) then
        transverse_parity = MOD(SUM(coord) - coord(mu), 2)
        colour = 2 * MOD(coord(mu), 2) + transverse_parity
     else
        colour = MOD(SUM(coord), 2)
     end if
   end function site_colour

   pure function site_linear_index(coord, dims) result(idx)
     !<
     !< Convert a 4D lattice coordinate to a unique 1-based linear site index.
     !<
     !< coord = [it, ix, iy, iz]
     !< dims  = [nt, nx, ny, nz]
     !<
     !< This linear index is used as counter(1) for the Philox stream so each
     !< site has its own independent random-number namespace.
     !<
     integer, intent(IN), dimension(4) :: coord
     !< site coordinate [it, ix, iy, iz]
     integer, intent(IN), dimension(4) :: dims
     !< total size of lattice [nt, nx, ny, nz]
     integer(kind=C64) :: idx
     !< linearised counter
      idx = INT((((coord(1) - 1) * dims(2) + (coord(2) - 1)) * dims(3) + &
                 (coord(3) - 1)) * dims(4) + (coord(4) - 1), C64) + 1_C64
   end function site_linear_index

end module FLUE_philox_helpers

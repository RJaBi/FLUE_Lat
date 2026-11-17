module FLUE
  !< Top level module exposing all external (i.e. app or python interface)
  !< subroutines/constants/functions
   use FLUE_constants, only: WP, WC, PI, SP, C_INT
   use FLUE_CSSM_bin, only: ReadGaugeTransformation_cola, ReadGaugeField_CSSM
   use FLUE_gluonProp, only: scalarGluonProp, calc_mom_space_scalarD
   use FLUE_gpManip, only: Q_Average, cone_cut
   use FLUE_heatbath, only: updateLinks, build_colour_sites
   use FLUE_ILDG_bin, only: ReadGaugeField_ILDG, writeGaugeField_ILDG
   use FLUE_ILDG_lime, only: writeGaugeField_ILDG_lime, writeGaugeField_ILDG_lime
   use FLUE_limeOut, only: writeLimeHeader, writeLimePadding
   use FLUE_matrixConstants, only: Ident3x3, Ident2x2, sigma1, sigma2, sigma3
   use FLUE_mom, only: get_qhat
   use FLUE_NERSCFileIO, only: WriteGaugeField_NERSC, readGaugeField_NERSC, header_NERSC, writeGaugeField_NERSC_FullHeader
   use FLUE_openQCDFileIO_SA, only: ReadGaugeField_OpenQCD, writeGaugeField_OpenQCD
   use FLUE_stoutSmearing, only: StoutSmearLinks
   use FLUE_SU2_CSSM, only: writeGaugeField_SU2_CSSM
   use FLUE_SU2_heatbath, only: SU2_updateLinks, constructXMatrix
   use FLUE_SU2_HKLS, only: ReadGaugeField_HKLS, WriteGaugeField_HKLS
   use FLUE_SU2_NRQ2CD, only: writeGaugeField_NRQ2CD, readGaugeField_NRQ2CD
   use FLUE_SU2_random, only: constructSU2Matrix
   use FLUE_SU2_wloops, only: SU2_genPlaquette, SU2_genericPath
   use FLUE_SU3_random, only: constructSU3Matrix
   use FLUE_SU3MatrixOps, only: MultiplyMatMat, MultiplyMatDagMatDag, &
                                TraceMultMatMat, RealTraceMultMatMat, TracelessConjgSubtract, &
                                colourDecomp, RealTraceMat, FixSU3Matrix, &
                                orthogonalise_vectors, vector_product, ExpIQ
   use FLUE_version, only: writeCompiler, writeGit
   use FLUE_wloops, only: polyakov, genPlaquette, magnetic, genericPath, periodCoord
   implicit none(type, external)

   character(len=*), parameter :: version = "0.1.4"
   !< Version of the FLUE software
   !public
   private
   public :: WP, WC
   !< constants
   public :: writeCompiler, writeGit
   public :: ReadGaugeField_CSSM
   public :: writeGaugeField_openQCD, readGaugeField_openQCD
   public :: ReadGaugeField_ILDG, writeGaugeField_ILDG
   public :: writeGaugeField_ILDG_Lime
   public :: writeGaugeField_NERSC
   public :: genPlaquette, magnetic
   public :: stoutSmearLinks
   public :: ReadGaugefield_HKLS, writeGaugeField_HKLS
   public :: writeGaugeField_SU2_CSSM
   public :: writeGaugeField_NRQ2CD
   public :: SU2_genPlaquette
   public :: Ident2x2, Ident3x3
   public :: Su2_updateLinks
   public :: updateLinks, build_colour_sites


end module FLUE

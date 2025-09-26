module FLUE

  use FLUE_SU2_HKLS, only: ReadGaugeField_HKLS, WriteGaugeField_HKLS
  use FLUE_Su2_CSSM, only: writeGaugeField_SU2_CSSM
  use FLUE_SU2_wloops, only: SU2_genPlaquette
  use FLUE_ILDG_bin, only: ReadGaugeField_ILDG, writeGaugeField_ILDG
  use FLUE_CSSM_bin, only: ReadGaugeTransformation_cola, ReadGaugeField_CSSM
   use FLUE_openQCDFileIO_SA, only: ReadGaugeField_OpenQCD, writeGaugeField_OpenQCD
   use FLUE_gluonProp, only: scalarGluonProp, calc_mom_space_scalarD
   use FLUE_wloops, only: plaquette, polyakov, genPlaquette, magnetic, genericPath, periodCoord
   use FLUE_SU3MatrixOps, only: Ident, MultiplyMatMat, MultiplyMatDagMatDag, &
        TraceMultMatMat, RealTraceMultMatMat, TracelessConjgSubtract, &
        colourDecomp, RealTraceMat, FixSU3Matrix, &
        orthogonalise_vectors, vector_product
   use FLUE_gpManip, only: Q_Average, cone_cut
   use FLUE_mom, only: get_qhat
   use FLUE_jack, only: complement, jackknife_wp
   use FLUE_constants, only: WP, WC, PI, SP, C_INT
   use FLUE_version, only: writeCompiler, writeGit
   implicit none(external)

   character(len=*), parameter :: version = "0.1.1"
   public
   !public :: calc_mom_space_scalarD
   !public :: Ident, MultiplyMatMat, MultiplyMatdagMatdag, TraceMultMatMat, RealTraceMultMatMat, TraceLessConjgSubtract, colourDecomp
   !public :: Q_Average
   !public :: get_qhat
   !public :: WP, pi, WC
   !public :: writeCompiler, writeGit

end module FLUE

module FLUE
   use FLUE_ILDG_bin, only: ReadGaugeField_ILDG, ReadGaugeTransformation_cola
   use FLUE_openQCDFileIO_SA, only: ReadGaugeField_OpenQCD
   use FLUE_gluonProp, only: scalarGluonProp, calc_mom_space_scalarD
   use FLUE_wloops, only: plaquette, polyakov, genPlaquette, magnetic
   use FLUE_SU3MatrixOps, only: Ident, MultiplyMatMat, MultiplyMatDagMatDag, TraceMultMatMat, RealTraceMultMatMat, &
        & TracelessConjgSubtract, colourDecomp, RealTraceMat
   use FLUE_gpManip, only: Q_Average, cone_cut
   use FLUE_mom, only: get_qhat
   use FLUE_constants, only: WP, WC, PI
   use FLUE_version, only: writeCompiler, writeGit
   implicit none(external)

   character(len=*), parameter :: version = "0.1.0"
   public
   !public :: calc_mom_space_scalarD
   !public :: Ident, MultiplyMatMat, MultiplyMatdagMatdag, TraceMultMatMat, RealTraceMultMatMat, TraceLessConjgSubtract, colourDecomp
   !public :: Q_Average
   !public :: get_qhat
   !public :: WP, pi, WC
   !public :: writeCompiler, writeGit

end module FLUE

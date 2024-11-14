module FLUE
  use FLUE_openQCDFileIO_SA
  use FLUE_gluonProp
  use FLUE_SU3MatrixOps
  use FLUE_gpManip
  use FLUE_mom
  use FLUE_constants
  use FLUE_version

  character(len=*), parameter :: version="0.1.0"
  !private
  !public :: calc_mom_space_scalarD
  !public :: Ident, MultiplyMatMat, MultiplyMatdagMatdag, TraceMultMatMat, RealTraceMultMatMat, TraceLessConjgSubtract, colourDecomp
  !public :: Q_Average
  !public :: get_qhat
  !public :: WP, pi, WC
  !public :: writeCompiler, writeGit
  
end module FLUE 

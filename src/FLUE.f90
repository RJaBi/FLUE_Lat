module FLUE
  use FLUE_gpManip
  use FLUE_mom
  use FLUE_constants
  use FLUE_version

  character(len=*), parameter :: version="0.1.0"
  !private
  !public :: version
  !public :: Q_Average
  !public :: get_qhat
  !public :: WP, pi, WC
  !public :: writeCompiler, writeGit
  
end module FLUE 

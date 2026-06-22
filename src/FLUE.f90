MODULE FLUE
   USE FLUE_constants, ONLY: WP, WC, PI, SP, C_INT
   USE FLUE_CSSM_bin, ONLY: ReadGaugeTransformation_cola, ReadGaugeField_CSSM
   USE FLUE_gluonProp, ONLY: scalarGluonProp, calc_mom_space_scalarD
   USE FLUE_gpManip, ONLY: Q_Average, cone_cut
   USE FLUE_heatbath, ONLY: updateLinks, build_colour_sites
   USE FLUE_ILDG_bin, ONLY: ReadGaugeField_ILDG, writeGaugeField_ILDG
   USE FLUE_matrixConstants, ONLY: Ident3x3, Ident2x2, sigma1, sigma2, sigma3
   USE FLUE_mom, ONLY: get_qhat
   USE FLUE_openQCDFileIO_SA, ONLY: ReadGaugeField_OpenQCD, writeGaugeField_OpenQCD
   USE FLUE_stoutSmearing, ONLY: StoutSmearLinks
   USE FLUE_SU2_CSSM, ONLY: writeGaugeField_SU2_CSSM
   USE FLUE_SU2_heatbath, ONLY: SU2_updateLinks, constructXMatrix
   USE FLUE_SU2_HKLS, ONLY: ReadGaugeField_HKLS, WriteGaugeField_HKLS
   USE FLUE_SU2_NRQ2CD, ONLY: writeGaugeField_NRQ2CD, readGaugeField_NRQ2CD
   USE FLUE_SU2_random, ONLY: constructSU2Matrix
   USE FLUE_SU2_wloops, ONLY: SU2_genPlaquette, SU2_genericPath
   USE FLUE_SU3_random, ONLY: constructSU3Matrix
   USE FLUE_SU3MatrixOps, ONLY: MultiplyMatMat, MultiplyMatDagMatDag, &
                                TraceMultMatMat, RealTraceMultMatMat, TracelessConjgSubtract, &
                                colourDecomp, RealTraceMat, FixSU3Matrix, &
                                orthogonalise_vectors, vector_product, ExpIQ
   USE FLUE_version, ONLY: writeCompiler, writeGit
   USE FLUE_wloops, ONLY: polyakov, genPlaquette, magnetic, genericPath, periodCoord
   IMPLICIT NONE(TYPE, EXTERNAL)

   CHARACTER(len=*), PARAMETER :: version = "0.1.1"
   PUBLIC
   !public :: calc_mom_space_scalarD
   !public :: Ident, MultiplyMatMat, MultiplyMatdagMatdag
   !public :: TraceMultMatMat, RealTraceMultMatMat, TraceLessConjgSubtract, colourDecomp
   !public :: Q_Average
   !public :: get_qhat
   !public :: WP, pi, WC
   !public :: writeCompiler, writeGit

END MODULE FLUE

########################################################################
#  This is the makefile for the eigenvalue test program from LAPACK.
#  The test files are organized as follows:
#
#     AEIGTST -- Auxiliary test routines used in all precisions
#     SCIGTST -- Auxiliary test routines used in single precision
#     DZIGTST -- Auxiliary test routines used in double precision
#     SEIGTST -- Single precision real test routines
#     CEIGTST -- Single precision complex test routines
#     DEIGTST -- Double precision real test routines
#     ZEIGTST -- Double precision complex test routines
#
########################################################################

set(AEIGTST
  alahdg.f alasum.f alasvm.f alareq.f
  ilaenv.f xerbla.f xlaenv.f chkxer.f)

set(SCIGTST
  slafts.f slahd2.f slasum.f slatb9.f sstech.f sstect.f
  ssvdch.f ssvdct.f ssxt1.f)

set(SEIGTST
  schkee.F
  sbdt01.f sbdt02.f sbdt03.f sbdt04.f sbdt05.f
  schkbb.f schkbd.f schkbk.f schkbl.f schkec.f
  schkgg.f schkgk.f schkgl.f schkhs.f schksb.f schkst.f schkst2stg.f schksb2stg.f
  sckcsd.f sckglm.f sckgqr.f sckgsv.f scklse.f scsdts.f
  sdrges.f sdrgev.f sdrges3.f sdrgev3.f sdrgsx.f sdrgvx.f
  sdrvbd.f sdrvmr3.f sdrves.f sdrvev.f sdrvsg.f sdrvsg2stg.f
  sdrvst.f sdrvst2stg.f sdrvsx.f sdrvvx.f
  serrbd.f serrec.f serred.f serrgg.f serrhs.f serrst.f
  sget02.f sget10.f sget22.f sget23.f sget24.f sget31.f
  sget32.f sget33.f sget34.f sget35.f sget36.f
  sget37.f sget38.f sget39.f sget40.f sget51.f sget52.f sget53.f
  sget54.f sglmts.f sgqrts.f sgrqts.f sgsvts3.f
  shst01.f slarfy.f slarhs.f slatm4.f slctes.f slctsx.f slsets.f sort01.f
  sort03.f ssbt21.f ssgt01.f sslect.f sspt21.f sstt21.f
  sstt22.f ssyl01.f ssyt21.f ssyt22.f)

set(SDMDEIGTST schkdmd.f90)

set(CEIGTST
  cchkee.F
  cbdt01.f cbdt02.f cbdt03.f cbdt05.f
  cchkbb.f cchkbd.f cchkbk.f cchkbl.f cchkec.f
  cchkgg.f cchkgk.f cchkgl.f cchkhb.f cchkhs.f cchkst.f cchkst2stg.f cchkhb2stg.f
  cckcsd.f cckglm.f cckgqr.f cckgsv.f ccklse.f ccsdts.f
  cdrges.f cdrgev.f cdrges3.f cdrgev3.f cdrgsx.f cdrgvx.f
  cdrvbd.f cdrvmr3.f cdrves.f cdrvev.f cdrvsg.f cdrvsg2stg.f
  cdrvst.f cdrvst2stg.f cdrvsx.f cdrvvx.f
  cerrbd.f cerrec.f cerred.f cerrgg.f cerrhs.f cerrst.f
  cget02.f cget10.f cget22.f cget23.f cget24.f
  cget35.f cget36.f cget37.f cget38.f cget51.f cget52.f
  cget54.f cglmts.f cgqrts.f cgrqts.f cgsvts3.f
  chbt21.f chet21.f chet22.f chpt21.f chst01.f
  clarfy.f clarhs.f clatm4.f clctes.f clctsx.f clsets.f csbmv.f
  csgt01.f cslect.f csyl01.f
  cstt21.f cstt22.f cunt01.f cunt03.f)

set(CDMDEIGTST cchkdmd.f90)

set(DZIGTST
  dlafts.f dlahd2.f dlasum.f dlatb9.f dstech.f dstect.f
  dsvdch.f dsvdct.f dsxt1.f)

set(DEIGTST
  dchkee.F
  dbdt01.f dbdt02.f dbdt03.f dbdt04.f dbdt05.f
  dchkbb.f dchkbd.f dchkbk.f dchkbl.f dchkec.f
  dchkgg.f dchkgk.f dchkgl.f dchkhs.f dchksb.f dchkst.f dchkst2stg.f dchksb2stg.f
  dckcsd.f dckglm.f dckgqr.f dckgsv.f dcklse.f dcsdts.f
  ddrges.f ddrgev.f ddrges3.f ddrgev3.f ddrgsx.f ddrgvx.f
  ddrvbd.f ddrvmr3.f ddrves.f ddrvev.f ddrvsg.f ddrvsg2stg.f
  ddrvst.f ddrvst2stg.f ddrvsx.f ddrvvx.f
  derrbd.f derrec.f derred.f derrgg.f derrhs.f derrst.f
  dget02.f dget10.f dget22.f dget23.f dget24.f dget31.f
  dget32.f dget33.f dget34.f dget35.f dget36.f
  dget37.f dget38.f dget39.f dget40.f dget51.f dget52.f dget53.f
  dget54.f dglmts.f dgqrts.f dgrqts.f dgsvts3.f
  dhst01.f dlarfy.f dlarhs.f dlatm4.f dlctes.f dlctsx.f dlsets.f dort01.f
  dort03.f dsbt21.f dsgt01.f dslect.f dspt21.f dstt21.f
  dstt22.f dsyl01.f dsyt21.f dsyt22.f)

set(DDMDEIGTST dchkdmd.f90)

set(ZEIGTST
  zchkee.F
  zbdt01.f zbdt02.f zbdt03.f zbdt05.f
  zchkbb.f zchkbd.f zchkbk.f zchkbl.f zchkec.f
  zchkgg.f zchkgk.f zchkgl.f zchkhb.f zchkhs.f zchkst.f zchkst2stg.f zchkhb2stg.f
  zckcsd.f zckglm.f zckgqr.f zckgsv.f zcklse.f zcsdts.f
  zdrges.f zdrgev.f zdrges3.f zdrgev3.f zdrgsx.f zdrgvx.f
  zdrvbd.f zdrvmr3.f zdrves.f zdrvev.f zdrvsg.f zdrvsg2stg.f
  zdrvst.f zdrvst2stg.f zdrvsx.f zdrvvx.f
  zerrbd.f zerrec.f zerred.f zerrgg.f zerrhs.f zerrst.f
  zget02.f zget10.f zget22.f zget23.f zget24.f
  zget35.f zget36.f zget37.f zget38.f zget51.f zget52.f
  zget54.f zglmts.f zgqrts.f zgrqts.f zgsvts3.f
  zhbt21.f zhet21.f zhet22.f zhpt21.f zhst01.f
  zlarfy.f zlarhs.f zlatm4.f zlctes.f zlctsx.f zlsets.f zsbmv.f
  zsgt01.f zslect.f zsyl01.f
  zstt21.f zstt22.f zunt01.f zunt03.f)

set(ZDMDEIGTST zchkdmd.f90)

macro(add_eig_executable name)
  add_executable(${name} ${ARGN})
  target_link_libraries(${name} ${LIBNAMEPREFIX}openblas${LIBNAMESUFFIX}${SUFFIX64_UNDERSCORE})
endmacro()

if(BUILD_SINGLE)
add_lapack_test(stest.out stest.in xlintsts)
#
# ======== SINGLE RFP LIN TESTS ========================
add_lapack_test(stest_rfp.out stest_rfp.in xlintstrfs)
#
#
# ======== SINGLE EIG TESTS ===========================
#
add_lapack_test(snep.out nep.in xeigtsts)
add_lapack_test(ssep.out sep.in xeigtsts)
add_lapack_test(sse2.out se2.in xeigtsts)
add_lapack_test(ssvd.out svd.in xeigtsts)
add_lapack_test(sec.out sec.in xeigtsts)
add_lapack_test(sed.out sed.in xeigtsts)
add_lapack_test(sgg.out sgg.in xeigtsts)
add_lapack_test(sgd.out sgd.in xeigtsts)
add_lapack_test(ssb.out ssb.in xeigtsts)
add_lapack_test(ssg.out ssg.in xeigtsts)
add_lapack_test(sbal.out sbal.in xeigtsts)
add_lapack_test(sbak.out sbak.in xeigtsts)
add_lapack_test(sgbal.out sgbal.in xeigtsts)
add_lapack_test(sgbak.out sgbak.in xeigtsts)
add_lapack_test(sbb.out sbb.in xeigtsts)
add_lapack_test(sglm.out glm.in xeigtsts)
add_lapack_test(sgqr.out gqr.in xeigtsts)
add_lapack_test(sgsv.out gsv.in xeigtsts)
add_lapack_test(scsd.out csd.in xeigtsts)
add_lapack_test(slse.out lse.in xeigtsts)
#
# ======== SINGLE DMD EIG TESTS ===========================
add_lapack_test(sdmd.out sdmd.in xdmdeigtsts)
add_test(NAME LAPACK-test_sgges3_scaling COMMAND ${CMAKE_CURRENT_BINARY_DIR}/EIG/test_sgges3_scaling)
add_test(NAME LAPACK-test_sgesvj_emptsw COMMAND ${CMAKE_CURRENT_BINARY_DIR}/EIG/test_sgesvj_emptsw)
endif()

if(BUILD_DOUBLE)
#
# ======== DOUBLE LIN TESTS ===========================
add_lapack_test(dtest.out dtest.in xlintstd)
#
# ======== DOUBLE RFP LIN TESTS ========================
add_lapack_test(dtest_rfp.out dtest_rfp.in xlintstrfd)
#
# ======== DOUBLE EIG TESTS ===========================
add_lapack_test(dnep.out nep.in xeigtstd)
add_lapack_test(dsep.out sep.in xeigtstd)
add_lapack_test(dse2.out se2.in xeigtstd)
add_lapack_test(dsvd.out svd.in xeigtstd)
add_lapack_test(dec.out dec.in xeigtstd)
add_lapack_test(ded.out ded.in xeigtstd)
add_lapack_test(dgg.out dgg.in xeigtstd)
add_lapack_test(dgd.out dgd.in xeigtstd)
add_lapack_test(dsb.out dsb.in xeigtstd)
add_lapack_test(dsg.out dsg.in xeigtstd)
add_lapack_test(dbal.out dbal.in xeigtstd)
add_lapack_test(dbak.out dbak.in xeigtstd)
add_lapack_test(dgbal.out dgbal.in xeigtstd)
add_lapack_test(dgbak.out dgbak.in xeigtstd)
add_lapack_test(dbb.out dbb.in xeigtstd)
add_lapack_test(dglm.out glm.in xeigtstd)
add_lapack_test(dgqr.out gqr.in xeigtstd)
add_lapack_test(dgsv.out gsv.in xeigtstd)
add_lapack_test(dcsd.out csd.in xeigtstd)
add_lapack_test(dlse.out lse.in xeigtstd)
#
# ======== DOUBLE DMD EIG TESTS ===========================
add_lapack_test(ddmd.out ddmd.in xdmdeigtstd)
add_test(NAME LAPACK-test_dgges3_scaling COMMAND ${CMAKE_CURRENT_BINARY_DIR}/EIG/test_dgges3_scaling)
add_test(NAME LAPACK-test_dgesvj_emptsw COMMAND ${CMAKE_CURRENT_BINARY_DIR}/EIG/test_dgesvj_emptsw)
endif()

if(BUILD_COMPLEX)
add_lapack_test(ctest.out ctest.in xlintstc)
#
# ======== COMPLEX RFP LIN TESTS ========================
add_lapack_test(ctest_rfp.out ctest_rfp.in xlintstrfc)
#
# ======== COMPLEX EIG TESTS ===========================
add_lapack_test(cnep.out nep.in xeigtstc)
add_lapack_test(csep.out sep.in xeigtstc)
add_lapack_test(cse2.out se2.in xeigtstc)
add_lapack_test(csvd.out svd.in xeigtstc)
add_lapack_test(cec.out cec.in xeigtstc)
add_lapack_test(ced.out ced.in xeigtstc)
add_lapack_test(cgg.out cgg.in xeigtstc)
add_lapack_test(cgd.out cgd.in xeigtstc)
add_lapack_test(csb.out csb.in xeigtstc)
add_lapack_test(csg.out csg.in xeigtstc)
add_lapack_test(cbal.out cbal.in xeigtstc)
add_lapack_test(cbak.out cbak.in xeigtstc)
add_lapack_test(cgbal.out cgbal.in xeigtstc)
add_lapack_test(cgbak.out cgbak.in xeigtstc)
add_lapack_test(cbb.out cbb.in xeigtstc)
add_lapack_test(cglm.out glm.in xeigtstc)
add_lapack_test(cgqr.out gqr.in xeigtstc)
add_lapack_test(cgsv.out gsv.in xeigtstc)
add_lapack_test(ccsd.out csd.in xeigtstc)
add_lapack_test(clse.out lse.in xeigtstc)
#
# ======== COMPLEX DMD EIG TESTS ===========================
add_lapack_test(cdmd.out cdmd.in xdmdeigtstc)
add_test(NAME LAPACK-test_cgges3_scaling COMMAND ${CMAKE_CURRENT_BINARY_DIR}/EIG/test_cgges3_scaling)
add_test(NAME LAPACK-test_cgesvj_emptsw COMMAND ${CMAKE_CURRENT_BINARY_DIR}/EIG/test_cgesvj_emptsw)
endif()

if(BUILD_COMPLEX16)
#
# ======== COMPLEX16 LIN TESTS ========================
add_lapack_test(ztest.out ztest.in xlintstz)
#
# ======== COMPLEX16 RFP LIN TESTS ========================
add_lapack_test(ztest_rfp.out ztest_rfp.in xlintstrfz)
#
# ======== COMPLEX16 EIG TESTS ===========================
add_lapack_test(znep.out nep.in xeigtstz)
add_lapack_test(zsep.out sep.in xeigtstz)
add_lapack_test(zse2.out se2.in xeigtstz)
add_lapack_test(zsvd.out svd.in xeigtstz)
add_lapack_test(zec.out zec.in xeigtstz)
add_lapack_test(zed.out zed.in xeigtstz)
add_lapack_test(zgg.out zgg.in xeigtstz)
add_lapack_test(zgd.out zgd.in xeigtstz)
add_lapack_test(zsb.out zsb.in xeigtstz)
add_lapack_test(zsg.out zsg.in xeigtstz)
add_lapack_test(zbal.out zbal.in xeigtstz)
add_lapack_test(zbak.out zbak.in xeigtstz)
add_lapack_test(zgbal.out zgbal.in xeigtstz)
add_lapack_test(zgbak.out zgbak.in xeigtstz)
add_lapack_test(zbb.out zbb.in xeigtstz)
add_lapack_test(zglm.out glm.in xeigtstz)
add_lapack_test(zgqr.out gqr.in xeigtstz)
add_lapack_test(zgsv.out gsv.in xeigtstz)
add_lapack_test(zcsd.out csd.in xeigtstz)
add_lapack_test(zlse.out lse.in xeigtstz)
#
# ======== COMPLEX16 DMD EIG TESTS ===========================
add_lapack_test(zdmd.out zdmd.in xdmdeigtstz)
add_test(NAME LAPACK-test_zgges3_scaling COMMAND ${CMAKE_CURRENT_BINARY_DIR}/EIG/test_zgges3_scaling)
add_test(NAME LAPACK-test_zgesvj_emptsw COMMAND ${CMAKE_CURRENT_BINARY_DIR}/EIG/test_zgesvj_emptsw)
endif()


if(BUILD_DOUBLE)
add_eig_executable(xeigtstd ${DEIGTST} ${DZIGTST} ${AEIGTST})
add_eig_executable(xdmdeigtstd ${DDMDEIGTST} ${AEIGTST})
add_eig_executable(test_dgges3_scaling test_dgges3_scaling.f90)
add_eig_executable(test_dgesvj_emptsw test_sgesvj_emptsw.f90)
endif()

if(BUILD_COMPLEX16)
add_eig_executable(xeigtstz ${ZEIGTST} ${DZIGTST} ${AEIGTST})
add_eig_executable(xdmdeigtstz ${ZDMDEIGTST} ${AEIGTST})
#add_eig_executable(test_wq_zrwork test_wq_rwork.f)
add_eig_executable(test_zgges3_scaling test_zgges3_scaling.f90)
add_eig_executable(test_zgesvj_emptsw test_sgesvj_emptsw.f90)
endif()

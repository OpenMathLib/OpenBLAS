# Sources for compiling lapack-netlib. Can't use CMakeLists.txt because lapack-netlib already has its own cmake files.

set(ALLAUX
  ilaenv.f ieeeck.f lsamen.f xerbla_array.f iparmq.f	
  ilaprec.f ilatrans.f ilauplo.f iladiag.f chla_transtype.f 
  ../INSTALL/ilaver.f ../INSTALL/slamch.f
)

set(DZLAUX
   dbdsdc.f 
   dbdsqr.f ddisna.f dlabad.f dlacpy.f dladiv.f dlae2.f  dlaebz.f 
   dlaed0.f dlaed1.f dlaed2.f dlaed3.f dlaed4.f dlaed5.f dlaed6.f 
   dlaed7.f dlaed8.f dlaed9.f dlaeda.f dlaev2.f dlagtf.f 
   dlagts.f dlamrg.f dlanst.f 
   dlapy2.f dlapy3.f dlarnv.f 
   dlarra.f dlarrb.f dlarrc.f dlarrd.f dlarre.f dlarrf.f dlarrj.f 
   dlarrk.f dlarrr.f dlaneg.f 
   dlartg.f dlaruv.f dlas2.f  dlascl.f 
   dlasd0.f dlasd1.f dlasd2.f dlasd3.f dlasd4.f dlasd5.f dlasd6.f 
   dlasd7.f dlasd8.f dlasda.f dlasdq.f dlasdt.f 
   dlaset.f dlasq1.f dlasq2.f dlasq3.f dlasq4.f dlasq5.f dlasq6.f 
   dlasr.f  dlasrt.f dlassq.f dlasv2.f dpttrf.f dstebz.f dstedc.f 
   dsteqr.f dsterf.f dlaisnan.f disnan.f 
   dlartgp.f dlartgs.f 
   ../INSTALL/dlamch.f ../INSTALL/dsecnd_${TIMER}.f
)

set(DLASRC
   dgbbrd.f dgbcon.f dgbequ.f dgbrfs.f dgbsv.f  
   dgbsvx.f dgbtf2.f dgbtrf.f dgbtrs.f dgebak.f dgebal.f dgebd2.f 
   dgebrd.f dgecon.f dgeequ.f dgees.f  dgeesx.f dgeev.f  dgeevx.f 
   dgegs.f  dgegv.f  dgehd2.f dgehrd.f dgelq2.f dgelqf.f 
   dgels.f  dgelsd.f dgelss.f dgelsx.f dgelsy.f dgeql2.f dgeqlf.f 
   dgeqp3.f dgeqpf.f dgeqr2.f dgeqr2p.f dgeqrf.f dgeqrfp.f dgerfs.f 
   dgerq2.f dgerqf.f dgesc2.f dgesdd.f  dgesvd.f dgesvx.f  
   dgetc2.f dgetri.f 
   dggbak.f dggbal.f dgges.f  dggesx.f dggev.f  dggevx.f 
   dggglm.f dgghrd.f dgglse.f dggqrf.f 
   dggrqf.f dggsvd.f dggsvp.f dgtcon.f dgtrfs.f dgtsv.f  
   dgtsvx.f dgttrf.f dgttrs.f dgtts2.f dhgeqz.f 
   dhsein.f dhseqr.f dlabrd.f dlacon.f dlacn2.f 
   dlaein.f dlaexc.f dlag2.f  dlags2.f dlagtm.f dlagv2.f dlahqr.f 
   dlahrd.f dlahr2.f dlaic1.f dlaln2.f dlals0.f dlalsa.f dlalsd.f 
   dlangb.f dlange.f dlangt.f dlanhs.f dlansb.f dlansp.f 
   dlansy.f dlantb.f dlantp.f dlantr.f dlanv2.f 
   dlapll.f dlapmt.f 
   dlaqgb.f dlaqge.f dlaqp2.f dlaqps.f dlaqsb.f dlaqsp.f dlaqsy.f 
   dlaqr0.f dlaqr1.f dlaqr2.f dlaqr3.f dlaqr4.f dlaqr5.f 
   dlaqtr.f dlar1v.f dlar2v.f iladlr.f iladlc.f 
   dlarf.f  dlarfb.f dlarfg.f dlarfgp.f dlarft.f dlarfx.f 
   dlargv.f dlarrv.f dlartv.f  
   dlarz.f  dlarzb.f dlarzt.f dlasy2.f dlasyf.f dlasyf_rook.f 
   dlatbs.f dlatdf.f dlatps.f dlatrd.f dlatrs.f dlatrz.f dlatzm.f 
   dopgtr.f dopmtr.f dorg2l.f dorg2r.f 
   dorgbr.f dorghr.f dorgl2.f dorglq.f dorgql.f dorgqr.f dorgr2.f 
   dorgrq.f dorgtr.f dorm2l.f dorm2r.f 
   dormbr.f dormhr.f dorml2.f dormlq.f dormql.f dormqr.f dormr2.f 
   dormr3.f dormrq.f dormrz.f dormtr.f dpbcon.f dpbequ.f dpbrfs.f 
   dpbstf.f dpbsv.f  dpbsvx.f 
   dpbtf2.f dpbtrf.f dpbtrs.f dpocon.f dpoequ.f dporfs.f dposv.f  
   dposvx.f dpotrs.f dpstrf.f dpstf2.f 
   dppcon.f dppequ.f 
   dpprfs.f dppsv.f  dppsvx.f dpptrf.f dpptri.f dpptrs.f dptcon.f 
   dpteqr.f dptrfs.f dptsv.f  dptsvx.f dpttrs.f dptts2.f drscl.f  
   dsbev.f  dsbevd.f dsbevx.f dsbgst.f dsbgv.f  dsbgvd.f dsbgvx.f 
   dsbtrd.f  dspcon.f dspev.f  dspevd.f dspevx.f dspgst.f 
   dspgv.f  dspgvd.f dspgvx.f dsprfs.f dspsv.f  dspsvx.f dsptrd.f 
   dsptrf.f dsptri.f dsptrs.f dstegr.f dstein.f dstev.f  dstevd.f dstevr.f 
   dstevx.f 
   dsycon.f dsyev.f  dsyevd.f dsyevr.f 
   dsyevx.f dsygs2.f dsygst.f dsygv.f  dsygvd.f dsygvx.f dsyrfs.f 
   dsysv.f  dsysvx.f 
   dsytd2.f dsytf2.f dsytrd.f dsytrf.f dsytri.f dsytri2.f dsytri2x.f 
   dsyswapr.f dsytrs.f dsytrs2.f dsyconv.f 
   dsytf2_rook.f dsytrf_rook.f dsytrs_rook.f 
   dsytri_rook.f dsycon_rook.f dsysv_rook.f 
   dtbcon.f dtbrfs.f dtbtrs.f dtgevc.f dtgex2.f dtgexc.f dtgsen.f 
   dtgsja.f dtgsna.f dtgsy2.f dtgsyl.f dtpcon.f dtprfs.f dtptri.f 
   dtptrs.f 
   dtrcon.f dtrevc.f dtrexc.f dtrrfs.f dtrsen.f dtrsna.f dtrsyl.f 
   dtrtrs.f dtzrqf.f dtzrzf.f dstemr.f 
   dsgesv.f dsposv.f dlag2s.f slag2d.f dlat2s.f 
   dlansf.f dpftrf.f dpftri.f dpftrs.f dsfrk.f dtfsm.f dtftri.f dtfttp.f 
   dtfttr.f dtpttf.f dtpttr.f dtrttf.f dtrttp.f 
   dgejsv.f  dgesvj.f  dgsvj0.f  dgsvj1.f 
   dgeequb.f dsyequb.f dpoequb.f dgbequb.f 
   dbbcsd.f dlapmr.f dorbdb.f dorbdb1.f dorbdb2.f dorbdb3.f dorbdb4.f 
   dorbdb5.f dorbdb6.f dorcsd.f dorcsd2by1.f 
   dgeqrt.f dgeqrt2.f dgeqrt3.f dgemqrt.f 
   dtpqrt.f dtpqrt2.f dtpmqrt.f dtprfb.f dpotri.f
)

set(LA_REL_SRC ${ALLAUX} ${DZLAUX} ${DLASRC})

# add lapack-netlib folder to the sources
set(LA_SOURCES "")
foreach (LA_FILE ${LA_REL_SRC})
  list(APPEND LA_SOURCES "${NETLIB_LAPACK_DIR}/SRC/${LA_FILE}")
endforeach ()


#
#  Activate CPPFLAGS depending on calculation type..  Default is MD/PIMD 
#
#
# === Makefile for DLPOLY with TTM2 ===
#
# Last updated: May 15, 2018 
# 

# target specific variables
#
#include machine/cab-intel
#include machine/vulcan-xlf
#include machine/cori-nersc-mpich-intel
#include machine/urey-openmpi-gfortran
#include machine/urey-openmpi-intel
#include machine/comet-mvapich2-intel
#include machine/stampede-mvapich2-intel
#include machine/gordon-mvapich2-intel
#include machine/tscc-openmpi
#include machine/dod-armstrong-mpich-intel-cray
#include machine/dod-conrad-mpich-intel-cray
#include machine/dod-gordon-mpich-intel-cray
#include machine/dod-shepard-mpich-intel-cray
#include machine/dod-thunder-mpich-intel
#include machine/dod-topaz-mpich-intel
#include machine/ubuntu
include machine/sissa-hpc-intel

# common variables
CC  = gcc 
CPP = cpp -P -traditional
CPPFLAGS += -DMPI -DSTRESS  # -DDEBUG
CPPFLAGS += -DHEAT_CURRENT  # for energy and stress decomposistion and 
                            #calculation of heat current
#CPPFLAGS += -DDO_2D_IR          # for 2D IR calculation 
#CPPFLAGS += -DOPTORVIB -DOPT     # for geometry optimization and vibrational analysis - DOPT is MB-pol to account for ceiling
#CPPFLAGS += -DDO_RECALC_DIP   # calculates DIPMOL_CMD, DIPIND_CMD for a given frame
#CPPFLAGS += -DIPI             #  interface to i-PI 
# for PLUMED, activate two lines below and also, define the PLUMEDHOME variable 
#CPPFLAGS += -DPLUMED             #  interface to PLUMED 
#LDFLAGS += "${PLUMEDHOME}/lib/libplumed.so" -ldl


DEBUG = 

#
# --- definitions of object files ---
#

OBJ_ALL = module.o heatcurrent.o angfrc.o bndfrc.o cfgscan.o corshl.o coul0.o coul4.o \
	coul2.o coul3.o conscan.o dblstr.o dcell.o diffsn0.o \
	diffsn1.o dlpoly.o duni.o error.o ewald1.o ewald3.o \
	exclude.o exclude_atom.o fldscan.o exclude_link.o \
	extnfld.o fbpfrc.o fcap.o forces.o freeze.o gauss.o \
	gdsum.o getrec.o gimax.o gisum.o gstate.o images.o initcomms.o \
	intlist.o intstr.o invert.o invfrc.o jacobi.o lowcase.o lrcmetal.o \
	lrcorrect.o merge.o merge1.o merge4.o multiple.o  \
	multiple_nsq.o npt_b1.o npt_b3.o parset.o npt_h1.o npt_h3.o nve_1.o \
	nvt_b1.o nvt_e1.o nvt_h1.o parlst_nsq.o parlink.o parlst.o passcon.o \
	passpmf.o pmf_1.o pmf_shake.o primlst.o quench.o rdf0.o rdf1.o \
	rdshake_1.o result.o revive.o scdens.o shellsort.o shlfrc.o \
	shlmerge.o shlqnch.o shmove.o simdef.o splice.o static.o strip.o \
	strucopt.o sysdef.o sysgen.o systemp.o sysbook.o sysinit.o \
	tethfrc.o thbfrc.o timchk.o traject.o traject_u.o vertest.o vscaleg.o \
	warning.o xscale.o zden0.o zden1.o \
	rfield.o erfcc.o ewald1p.o iter.o coulomb_polar.o ewald2p.o ewald3p.o \
	ewald4p.o erfcgenp.o ewald1cp.o ewald2cp.o ewald3cp.o ewald4cp.o \
	thole.o dcp_0.o dcp_h0.o cenmas.o \
	geopt_pt.o engforce.o frprmn.o ran2.o \
	ttm2list.o gammq.o output_ttm2_dip.o \
	output_ttm2_com.o readrec.o fordebug.o ps_type_couple.o checkcom.o \
	output_ttm2_force.o qdetermine.o qdforce.o optimization.o \
	set_guess_dipole.o multibead.o exclude_ttm.o tang-toennies.o aspc.o

OBJ_SPME = bspcoe.o bspgen.o cpy_rtc.o ele_prd.o ewald_spme.o \
	scl_csum.o set_block.o spl_cexp.o spme_for.o dlpfft3.o

OBJ_HKE = hkgen.o hkewald1.o hkewald2.o hkewald3.o hkewald4.o cerfr.o

OBJ_NEU = coul0neu.o coul2neu.o coul3neu.o excludeneu.o forcesneu.o \
	multipleneu.o neutlst.o parneulst.o prneulst.o \
	parlinkneu.o rdf0neu.o 

OBJ_RIG = nptq_b1.o nptq_b2.o nptq_b3.o nptq_b4.o nptq_h1.o nptq_h2.o \
	nptq_h3.o nptq_h4.o nveq_1.o nveq_2.o nvtq_b1.o nvtq_b2.o \
	nvtq_h1.o nvtq_h2.o passquat.o qshake.o quatbook.o quatqnch.o

OBJ_RRR = denloc.o dihfrc.o erfcgen.o ewald2.o ewald4.o forgen.o \
	fortab.o metgen.o srfrce.o srfrceneu.o suttchen.o srfrce2.o

OBJ_4PT = denloc_4pt.o dihfrc_4pt.o erfcgen.o ewald2_4pt.o ewald4_4pt.o \
	forgen.o fortab.o metgen.o srfrce_4pt.o srfrceneu_4pt.o \
	suttchen_4pt.o

OBJ_RSQ = denloc_rsq.o dihfrc_rsq.o erfcgen_rsq.o ewald2_rsq.o ewald4_rsq.o \
	forgen_rsq.o fortab_rsq.o metgen_rsq.o srfrce_rsq.o srfrceneu_rsq.o\
	suttchen_rsq.o 

OBJ_EXT = crecv.o csend.o gsync.o nodedim.o pdsum.o assert.o fixme.o

TIMER = etime.o

#
# --- pick up necessary objects ---
#
# 3pt 
OBJS = $(OBJ_ALL) $(OBJ_RRR) $(OBJ_NEU) $(OBJ_RIG) $(OBJ_EXT) $(TIMER) \
	$(OBJ_SPME) $(OBJ_HKE) 

# 4pt
#OBJS = $(OBJ_ALL) $(OBJ_4PT) $(OBJ_NEU) $(OBJ_RIG) $(OBJ_EXT) $(TIMER) \
#	$(OBJ_SPME) $(OBJ_HKE) 

# r-squared
#OBJS = $(OBJ_ALL) $(OBJ_RSQ)  $(OBJ_NEU) $(OBJ_RIG) $(OBJ_EXT) $(TIMER) \
#	$(OBJ_SPME) $(OBJ_HKE)

CMD_DLPOLY_OBJS = \
angfrc.o bndfrc.o checkcom.o corshl.o dcell.o dcp_0.o dcp_h0.o \
dihfrc.o engforce.o error.o extnfld.o fbpfrc.o fcap.o forces.o \
forcesneu.o fordebug.o freeze.o gdsum.o gstate.o images.o invfrc.o \
lrcmetal.o multiple.o multipleneu.o parlink.o parlinkneu.o parlst.o \
parneulst.o parset.o qdetermine.o quatqnch.o result.o revive.o shlfrc.o \
shlqnch.o static.o strucopt.o sysbook.o sysdef.o \
sysinit.o tethfrc.o thbfrc.o timchk.o traject.o traject_u.o \
ttm2list.o vertest.o vscaleg.o xscale.o zden0.o module.o heatcurrent.o gsync.o \
getrec.o intstr.o dblstr.o gimax.o gisum.o pdsum.o merge.o \
ewald1.o ewald1p.o ewald1cp.o ewald2.o ewald2p.o ewald2cp.o \
ewald3.o ewald3p.o ewald3cp.o ewald4.o ewald4p.o ewald4cp.o \
coulomb_polar.o qdforce.o fldscan.o cfgscan.o \
hkewald1.o hkewald2.o hkewald3.o ewald_spme.o rfield.o \
coul0.o coul2.o coul3.o coul4.o lowcase.o strip.o erfcc.o bspcoe.o \
forgen.o fortab.o readrec.o spme_for.o erfcgen.o invert.o gauss.o \
iter.o rdf0.o hkgen.o warning.o duni.o coul0neu.o coul2neu.o coul3neu.o \
merge1.o dlpfft3.o srfrce.o srfrce2.o suttchen.o scdens.o neutlst.o \
rdf0neu.o srfrceneu.o primlst.o hkewald4.o prneulst.o shellsort.o \
conscan.o rdf1.o diffsn0.o diffsn1.o shlmerge.o quatbook.o zden1.o \
exclude.o intlist.o passquat.o excludeneu.o exclude_link.o metgen.o \
erfcgenp.o quench.o gammq.o bspgen.o set_block.o cpy_rtc.o lrcorrect.o \
ele_prd.o scl_csum.o merge4.o jacobi.o denloc.o passcon.o passpmf.o \
exclude_atom.o shmove.o splice.o ps_type_couple.o parlst_nsq.o \
multiple_nsq.o pmf_1.o rdshake_1.o pmf_shake.o set_guess_dipole.o \
output_ttm2_force.o output_ttm2_com.o output_ttm2_dip.o optimization.o \
qshake.o nve_1.o nveq_1.o nveq_2.o nptq_h3.o nptq_h4.o npt_b1.o npt_b3.o \
npt_h1.o nvt_h1.o nvt_b1.o nvt_e1.o nptq_b1.o nptq_b2.o nptq_b3.o nptq_b4.o \
nptq_h1.o nptq_h2.o nvtq_b1.o nvtq_b2.o nvtq_h1.o nvtq_h2.o npt_h3.o \
geopt_pt.o frprmn.o cenmas.o ran2.o exclude_ttm.o tang-toennies.o aspc.o \
langevin.o cv.o firstlast_field.o optimization_module.o lbfgs.o \
vibrational_module.o initcomms.o 

MBPOL_OBJS = mbpol.o mbpol-2b-poly.o poly-2b-v6x.o \
                     mbpol-3b-poly.o poly-3b-v2x.o

MBNRG_OBJS = mbnrg.o mbnrg-2b-h2o-f-poly.o \
                     mbnrg-2b-h2o-cl-poly.o \
                     mbnrg-2b-h2o-br-poly.o \
                     mbnrg-2b-h2o-i-poly.o \
                     mbnrg-2b-h2o-li-poly.o \
                     mbnrg-2b-h2o-na-poly.o \
                     mbnrg-2b-h2o-k-poly.o \
                     mbnrg-2b-h2o-rb-poly.o \
                     mbnrg-2b-h2o-cs-poly.o \
                     mbnrg-2b-f-f-poly.o \
                     mbnrg-2b-cl-cl-poly.o \
                     mbnrg-2b-br-br-poly.o \
                     mbnrg-2b-i-i-poly.o \
                     mbnrg-2b-li-li-poly.o \
                     mbnrg-2b-na-na-poly.o \
                     mbnrg-2b-k-k-poly.o \
                     mbnrg-2b-rb-rb-poly.o \
                     mbnrg-2b-cs-cs-poly.o \
                     mbnrg-2b-li-f-poly.o \
                     mbnrg-2b-li-cl-poly.o \
                     mbnrg-2b-li-br-poly.o \
                     mbnrg-2b-li-i-poly.o \
                     mbnrg-2b-na-f-poly.o \
                     mbnrg-2b-na-cl-poly.o \
                     mbnrg-2b-na-br-poly.o \
                     mbnrg-2b-na-i-poly.o \
                     mbnrg-2b-k-f-poly.o \
                     mbnrg-2b-k-cl-poly.o \
                     mbnrg-2b-k-br-poly.o \
                     mbnrg-2b-k-i-poly.o \
                     mbnrg-2b-rb-f-poly.o \
                     mbnrg-2b-rb-cl-poly.o \
                     mbnrg-2b-rb-br-poly.o \
                     mbnrg-2b-rb-i-poly.o \
                     mbnrg-2b-cs-f-poly.o \
                     mbnrg-2b-cs-cl-poly.o \
                     mbnrg-2b-cs-br-poly.o \
                     mbnrg-2b-cs-i-poly.o \
                     mbnrg-3b-h2o-h2o-f-poly.o \
                     mbnrg-3b-h2o-h2o-cl-poly.o \
                     mbnrg-3b-h2o-h2o-br-poly.o \
            	     mbnrg-3b-h2o-h2o-i-poly.o \
                     mbnrg-3b-h2o-h2o-li-poly.o \
                     mbnrg-3b-h2o-h2o-na-poly.o \
                     mbnrg-3b-h2o-h2o-k-poly.o \
                     mbnrg-3b-h2o-h2o-rb-poly.o \
                     mbnrg-3b-h2o-h2o-cs-poly.o \
                     poly-2b-h2o-ion-v1x.o \
                     poly-2b-A1-A1-v1x.o \
                     poly-2b-A1-B1-v1x.o \
                     poly-3b-h2o-ion-v1x.o 

#mbnrg-3b-h2o-ion-poly.o 

CAVITY_OBJS = centroid_make_cavity.o

CAVITY_OBJS = centroid_make_cavity.o

CMD_OBJS = $(CMD_DLPOLY_OBJS) simdef_cmd.o sysdef_cmd.o dlpoly_cmd.o \
sysgen_cmd.o systemp_cmd.o multibead.o assert.o centroid_module.o \
nose_hoover_module.o centroid_vertest.o centroid_spring_force.o \
centroid_write_config.o centroid_write_dipole.o centroid_traject_cmd.o \
centroid_read_config.o centroid_init.o centroid_matrix.o \
centroid_read_dipole.o fixme.o centroid_traject_pimd.o \
induced_error.o centroid_dcp.o centroid_2D_IR.o centroid_recalc_dip.o \
centroid_traject_bead.o $(CAVITY_OBJS) $(MBPOL_OBJS) $(MBNRG_OBJS) 

IPI_OBJS = $(CMD_OBJS) ipi.o sockets.o fsockets.o

#
# --- program link ---
#

$(PROGRAM): $(OBJS) 
# 	@echo "Linking $(PROGRAM) ..."
#	$(LD) $(OBJS) $(LDFLAGS) $(LDLIBS) -o $(PROGRAM) 
#	@echo "done"

PROGRAM_CMD = $(PROGRAM)_cmd

cmd: $(PROGRAM_CMD)

$(PROGRAM_CMD): $(CMD_OBJS)
	$(LD) $(CMD_OBJS) $(LDFLAGS) $(LDLIBS) -o $(PROGRAM_CMD)


PROGRAM_IPI = $(PROGRAM)_ipi

ipi: $(PROGRAM_IPI)

$(PROGRAM_IPI): $(IPI_OBJS)
	$(LD) $(IPI_OBJS) $(LDFLAGS) $(LDLIBS) -o $(PROGRAM_IPI)


PROGRAM_2DIR = $(PROGRAM)_2dir

2dir: $(PROGRAM) $(PROGRAM_2DIR)


$(PROGRAM_2DIR): $(CMD_OBJS)
	$(LD) $(CMD_OBJS) $(LDFLAGS) $(LDLIBS) -o $(PROGRAM_2DIR)


PROGRAM_MAPDIP = $(PROGRAM)_mapdip

mapdip: $(PROGRAM) $(PROGRAM_MAPDIP)


$(PROGRAM_MAPDIP): $(CMD_OBJS)
	$(LD) $(CMD_OBJS) $(LDFLAGS) $(LDLIBS) -o $(PROGRAM_MAPDIP)



#
# --- clean up ---
#

clean:
	rm -f $(OBJ_ALL) $(OBJ_RRR) $(OBJ_EXT) $(OBJ_NEU) \
        $(OBJ_RIG) $(TIMER) $(OBJ_SPME) $(OBJ_4PT) $(OBJ_RSQ) $(OBJ_HKE)\
        $(CMD_OBJS) $(IPI_OBJS)   _*.f _*.f90 mpif.h *.mod *~

clean-ipi:
	rm -f ipi.o sockets.o simdef_cmd.o centroid_init.o static.o dlpoly_cmd.o 

#
# --- dependencies ---
#

.SUFFIXES:

etime.o: etime.c
	$(CC) -c etime.c

sockets.o: sockets.c
	$(CC) -c sockets.c

include depend

include mbpol.mk

.PHONY: clean cmd ipi 2dir mapdip

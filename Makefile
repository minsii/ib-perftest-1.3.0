TESTS = write_bw_mt_nqp write_bw_mt_nctx_iter write_bw_mt_nqp_ncq_qpiter write_bw_mt_nqp_qpiter write_mr_mt_1qp
UTILS = clock_test
MPI_TEST=write_mr_mpi_np_1qp

all: ${RDMACM_TESTS} ${MCAST_TESTS} ${TESTS} ${UTILS}
CC=icc
MPICC=mpicc
MPI_CFLAGS += -D_GNU_SOURCE -O2 -Wall
CFLAGS += -D_GNU_SOURCE -O2 -openmp -D_OPENMP_ON -Wall
OPENMP_CFLAGS = -openmp $(CFLAGS)
CFLAGS += -I/opt/intel/mic/ofed/card/usr/include -L/opt/intel/mic/ofed/card/usr/lib64
BASIC_FILES = get_clock.c
EXTRA_FILES = perftest_resources.c perftest_communication.c perftest_parameters.c
MCAST_FILES = multicast_resources.c
BASIC_HEADERS = get_clock.h
EXTRA_HEADERS = perftest_resources.h perftest_communication.h perftest_parameters.h
MCAST_HEADERS = multicast_resources.h
#The following seems to help GNU make on some platforms
LOADLIBES += -libverbs -lrdmacm
LDFLAGS +=

#${MCAST_TESTS}: LOADLIBES += -libumad -lm
	
${RDMACM_TESTS}: %: %.c ${BASIC_FILES} ${BASIC_HEADERS}
	$(CC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) $< ${BASIC_FILES} $(LOADLIBES) $(LDLIBS) -o $@
#${MCAST_TESTS}: %: %.c ${BASIC_FILES} ${EXTRA_FILES} ${MCAST_FILES} ${BASIC_HEADERS} ${EXTRA_HEADERS} ${MCAST_HEADERS}
#	$(CC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) $< ${BASIC_FILES} ${EXTRA_FILES} ${MCAST_FILES} $(LOADLIBES) $(LDLIBS) -o ib_$@
${TESTS} ${UTILS}: %: %.c ${BASIC_FILES} ${EXTRA_FILES} ${BASIC_HEADERS} ${EXTRA_HEADERS}
	$(CC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) $< ${BASIC_FILES} ${EXTRA_FILES} $(LOADLIBES) $(LDLIBS) -o ib_$@
${MPI_TEST}: %: %.c ${BASIC_FILES} ${EXTRA_FILES} ${BASIC_HEADERS} ${EXTRA_HEADERS}
	$(MPICC) $(CPPFLAGS) $(MPI_CFLAGS) $(LDFLAGS) $< ${BASIC_FILES} ${EXTRA_FILES} $(LOADLIBES) $(LDLIBS) -o ib_$@

clean:
	$(foreach fname,${RDMACM_TESTS}, rm -f ${fname})
	$(foreach fname,${MCAST_TESTS}, rm -f ib_${fname})
	$(foreach fname,${TESTS} ${UTILS}, rm -f ib_${fname})
	$(foreach fname,${MPI_TEST}, rm -f ib_${fname})
.DELETE_ON_ERROR:
.PHONY: all clean

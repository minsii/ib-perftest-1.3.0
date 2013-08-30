#RDMACM_TESTS = rdma_lat rdma_bw
#MCAST_TESTS = send_bw send_lat
#TESTS = write_bw_postlist write_lat write_bw read_lat read_bw atomic_lat atomic_bw
#TESTS = write_lat write_bw
TESTS = write_bw_mt_nqp
UTILS = clock_test


all: ${RDMACM_TESTS} ${MCAST_TESTS} ${TESTS} ${UTILS}
CC=icc
#CFLAGS += -Wall -g -D_GNU_SOURCE -O2 -mmic
CFLAGS += -Wall -g -D_GNU_SOURCE -O2 -mmic -openmp -D_OPENMP_ON
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

clean:
	$(foreach fname,${RDMACM_TESTS}, rm -f ${fname})
	$(foreach fname,${MCAST_TESTS}, rm -f ib_${fname})
	$(foreach fname,${TESTS} ${UTILS}, rm -f ib_${fname})
.DELETE_ON_ERROR:
.PHONY: all clean

/*
 * Copyright (c) 2005 Topspin Communications.  All rights reserved.
 * Copyright (c) 2005 Mellanox Technologies Ltd.  All rights reserved.
 * Copyright (c) 2009 HNR Consulting.  All rights reserved.
 * Copyright (c) 2013 by Argonne National Laboratory.  All rights reserved.
 *
 * This software is available to you under a choice of one of two
 * licenses.  You may choose to be licensed under the terms of the GNU
 * General Public License (GPL) Version 2, available from the file
 * COPYING in the main directory of this source tree, or the
 * OpenIB.org BSD license below:
 *
 *     Redistribution and use in source and binary forms, with or
 *     without modification, are permitted provided that the following
 *     conditions are met:
 *
 *      - Redistributions of source code must retain the above
 *        copyright notice, this list of conditions and the following
 *        disclaimer.
 *
 *      - Redistributions in binary form must reproduce the above
 *        copyright notice, this list of conditions and the following
 *        disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
 * ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 * $Id$
 */


/**
 * OpenMP version
 * This experiment suite proves that multi-threading can improve BW for small messages.
 *
 * This test parallels ibv_post_send with shared CTX, separate QP, shared CQ.
 * Each thread issues ibv_post_send for different QP in parallel, and polls CQ in sequential.
 */

#if HAVE_CONFIG_H
#  include <config.h>
#endif /* HAVE_CONFIG_H */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <malloc.h>
#include <sys/errno.h>

#include "get_clock.h"
#include "perftest_parameters.h"
#include "perftest_resources.h"
#include "perftest_communication.h"


//#define DEBUG
#ifdef DEBUG
#define debug_printf(fmt...) {printf(fmt);fflush(stdout);}
#else
#define debug_printf(fmt...)
#endif

#define VERSION 2.6

cycles_t	*tposted;
cycles_t	*tcompleted;

/****************************************************************************** 
 *
 ******************************************************************************/
static int pp_connect_ctx(struct pingpong_context *ctx,int my_psn,
						  struct pingpong_dest *dest, 
						  struct perftest_parameters *user_parm, int qpindex)
{
	struct ibv_qp_attr attr;
	memset(&attr, 0, sizeof attr);

	attr.qp_state 		= IBV_QPS_RTR;
	attr.path_mtu       = user_parm->curr_mtu;
	attr.dest_qp_num 	= dest->qpn;
	attr.rq_psn 		= dest->psn;
	attr.ah_attr.dlid   = dest->lid;
	if (user_parm->connection_type==RC) {
		attr.max_dest_rd_atomic     = 1;
		attr.min_rnr_timer          = 12;
	}
	if (user_parm->gid_index<0) {
		attr.ah_attr.is_global  = 0;
		attr.ah_attr.sl         = user_parm->sl;
	} else {
		attr.ah_attr.is_global  = 1;
		attr.ah_attr.grh.dgid   = dest->gid;
		attr.ah_attr.grh.sgid_index = user_parm->gid_index;
		attr.ah_attr.grh.hop_limit = 1;
		attr.ah_attr.sl         = 0;
	}
	attr.ah_attr.src_path_bits = 0;
	attr.ah_attr.port_num   = user_parm->ib_port;
	if (user_parm->connection_type == RC) {
		if (ibv_modify_qp(ctx->qp[qpindex], &attr,
				  IBV_QP_STATE              |
				  IBV_QP_AV                 |
				  IBV_QP_PATH_MTU           |
				  IBV_QP_DEST_QPN           |
				  IBV_QP_RQ_PSN             |
				  IBV_QP_MIN_RNR_TIMER      |
				  IBV_QP_MAX_DEST_RD_ATOMIC)) {
			fprintf(stderr, "Failed to modify RC QP to RTR\n");
			return 1;
		}
		attr.timeout            = user_parm->qp_timeout;
		attr.retry_cnt          = 7;
		attr.rnr_retry          = 7;
	} else {
		if (ibv_modify_qp(ctx->qp[qpindex], &attr,
				  IBV_QP_STATE              |
				  IBV_QP_AV                 |
				  IBV_QP_PATH_MTU           |
				  IBV_QP_DEST_QPN           |
				  IBV_QP_RQ_PSN)) {
			fprintf(stderr, "Failed to modify UC QP to RTR\n");
			return 1;
		}

	}
	attr.qp_state 	    = IBV_QPS_RTS;
	attr.sq_psn 	    = my_psn;
	attr.max_rd_atomic  = 1;
	if (user_parm->connection_type == 0) {
		attr.max_rd_atomic  = 1;
		if (ibv_modify_qp(ctx->qp[qpindex], &attr,
				  IBV_QP_STATE              |
				  IBV_QP_SQ_PSN             |
				  IBV_QP_TIMEOUT            |
				  IBV_QP_RETRY_CNT          |
				  IBV_QP_RNR_RETRY          |
				  IBV_QP_MAX_QP_RD_ATOMIC)) {
			fprintf(stderr, "Failed to modify RC QP to RTS\n");
			return 1;
		}
	} else {
		if (ibv_modify_qp(ctx->qp[qpindex], &attr,
				  IBV_QP_STATE              |
				  IBV_QP_SQ_PSN)) {
			fprintf(stderr, "Failed to modify UC QP to RTS\n");
			return 1;
		}

	}
	return 0;
}

/****************************************************************************** 
 *
 ******************************************************************************/
static void print_report(struct perftest_parameters *user_param) {

	double cycles_to_units;
	unsigned long tsize;	/* Transferred size, in megabytes */
	int i, j;
	int opt_posted = 0, opt_completed = 0;
	cycles_t opt_delta;
	cycles_t t;
	int iters = user_param->iters;
	int nthreads;

#pragma omp parallel
	nthreads = omp_get_num_threads();

	opt_delta = tcompleted[opt_posted] - tposted[opt_completed];

//	for (j = 0; j < iters * user_param->num_of_qps; j+=iters) {
	    debug_printf("%lld %lld\n", tcompleted[(iters*user_param->num_of_qps) - 1], tposted[0]);
//	}

	if (user_param->noPeak == OFF && 0) {
		/* Find the peak bandwidth unless asked not to in command line*/
		for (i = 0; i < iters * user_param->num_of_qps; ++i)
		  for (j = i; j < iters * user_param->num_of_qps; ++j) {
		    t = (tcompleted[j] - tposted[i]) / (j - i + 1);
		    if (t < opt_delta) {
		      opt_delta  = t;
		      opt_posted = i;
		      opt_completed = j;
		    }
		  }
	}
	
	cycles_to_units = get_cpu_mhz(user_param->cpu_freq_f) * 1000000;

	tsize = user_param->duplex ? 2 : 1;
	tsize = tsize * user_param->size;
	printf("%d %lu %d %.2f\n", nthreads, (unsigned long)user_param->size, iters,
	       tsize*(iters - user_param->tx_depth)*user_param->num_of_qps*cycles_to_units/(tcompleted[(iters*user_param->num_of_qps) - 1] - tposted[0]) / 0x100000);
}

/****************************************************************************** 
 *
 ******************************************************************************/
inline void increase_rem_addr_ncnt(struct ibv_send_wr *wr, int size, int scnt,
        uint64_t prim_addr, VerbType verb) {

    if (verb == ATOMIC) {
        wr->wr.atomic.remote_addr = prim_addr + INC(size) * scnt % CYCLE_BUFFER;
    } else {
        wr->wr.rdma.remote_addr = prim_addr + INC(size) * scnt % CYCLE_BUFFER;
    }
}

/******************************************************************************
 *
 ******************************************************************************/
inline void increase_loc_addr_ncnt(struct ibv_sge *sg,int size,int rcnt,
                              uint64_t prim_addr,int server_is_ud) {

    sg->addr = prim_addr +  INC(size) * rcnt % CYCLE_BUFFER;
}

#undef DEF_WC_SIZE
#define DEF_WC_SIZE         (300)

/******************************************************************************
 *
 ******************************************************************************/
int run_iter(struct pingpong_context *ctx, 
			 struct perftest_parameters *user_param,
			 struct pingpong_dest *rem_dest)	{

    int                totscnt = 0;
	int 			   totccnt = 0;
	int                i       = 0;
    int                index,ne;
	int				   warmindex;
    struct ibv_wc 	   *wc       = NULL;
	struct ibv_sge     *sge_list = NULL;
    struct ibv_send_wr *wr       = NULL;
	uint64_t		   *my_addr  = NULL;
	uint64_t		   *rem_addr = NULL;
	int *iter_depth; /* available depth of each iteration */

	ALLOCATE(wr ,struct ibv_send_wr , user_param->num_of_qps);
	ALLOCATE(sge_list ,struct ibv_sge , user_param->num_of_qps);
	ALLOCATE(my_addr ,uint64_t ,user_param->num_of_qps);
	ALLOCATE(rem_addr ,uint64_t ,user_param->num_of_qps);
	ALLOCATE(wc ,struct ibv_wc , DEF_WC_SIZE);
	ALLOCATE(iter_depth ,int ,user_param->num_of_qps);

	// Each QP has its own wr and sge , that holds the qp addresses and attr.
	// We write in cycles on the buffer to exploid the "Nahalem" system.
#pragma omp parallel for
	for (index = 0 ; index < user_param->num_of_qps ; index++) {

		sge_list[index].addr   = (uintptr_t)ctx->buf + (index*BUFF_SIZE(ctx->size));
		sge_list[index].length = user_param->size;
		sge_list[index].lkey   = ctx->mr->lkey;

		wr[index].sg_list             = &sge_list[index]; 
		wr[index].num_sge             = MAX_SEND_SGE;
		wr[index].opcode	          = IBV_WR_RDMA_WRITE;
		wr[index].next                = NULL;
		wr[index].wr.rdma.remote_addr = rem_dest[index].vaddr;
		wr[index].wr.rdma.rkey        = rem_dest[index].rkey;
		wr[index].wr_id               = index + 1;
		wr[index].send_flags          = IBV_SEND_SIGNALED;

		if (user_param->size <= user_param->inline_size) 
			wr[index].send_flags |= IBV_SEND_INLINE;

		ctx->scnt[index] = 0;
		ctx->ccnt[index] = 0;
		my_addr[index]	 = sge_list[index].addr;
		rem_addr[index]  = wr[index].wr.rdma.remote_addr;

		iter_depth[index] = user_param->tx_depth;
	}
	
	int errs = 0;

#ifdef _OPENMP_ON
#pragma omp parallel for reduction(+: errs) shared(wr, rem_addr, my_addr) private(warmindex, index)
#endif
	for (index = 0; index < user_param->num_of_qps; index++) {
	    for (warmindex = 0; warmindex < user_param->tx_depth; warmindex++) {
            /* private per thread */
            struct ibv_send_wr _wr = wr[index];
            struct ibv_send_wr *bad_wr;
            _wr.send_flags |= IBV_SEND_SIGNALED;
            _wr.wr_id = index + 1;

            if(_wr.wr_id != index + 1)
                fprintf(stderr, "1)wr_id id wrong %ld != %d\n", _wr.wr_id, index + 1);

            if (ibv_post_send(ctx->qp[index], &_wr, &bad_wr)) {
                fprintf(stderr, "1)Couldn't post send: (depth:%d, qp:%d), "
                        "wr num_seg %d, addr 0x%lx, length %d, lkey 0x%x, flag 0x%x \n",
                        warmindex, index, _wr.num_sge, _wr.sg_list[0].addr,
                        _wr.sg_list[0].length, _wr.sg_list[0].lkey, _wr.send_flags);
                errs++;
            }

            if(_wr.wr_id != index + 1)
                fprintf(stderr, "2)wr_id id wrong %ld != %d\n", _wr.wr_id, index + 1);
        }
	}
    if(errs > 0){
        return 1;
    }

    // set qp full
    for (index = 0; index < user_param->num_of_qps; index++) {
        iter_depth[index] = 0;
    }
	totscnt = user_param->tx_depth * user_param->num_of_qps;

	// clear cq
    do {
        ne = ibv_poll_cq(ctx->send_cq, DEF_WC_SIZE, wc);
        if (ne > 0) {
            for (i = 0; i < ne; i++) {

                if (wc[i].status != IBV_WC_SUCCESS)
                    NOTIFY_COMP_ERROR_SEND(wc[i],totscnt,totccnt);

                iter_depth[(int)wc[i].wr_id - 1]++; // set empty
                totccnt++;
            }
        }
        else if(ne < 0){
            break;
        }
    } while (ne > 0 && totccnt < totscnt);
    if (ne < 0) {
        fprintf(stderr, "poll CQ failed %d\n", ne);
        return 1;
    }
    tcompleted[totccnt-1] = get_cycles();

    debug_printf("1-2) totscnt %d totccnt %d\n", totscnt, totccnt);
#ifdef DEBUG
    for (index = 0; index < user_param->num_of_qps; index++) {
        debug_printf("depth[%d] = %d \n", index, iter_depth[index]);
    }
#endif

    int main_iter = 2;
    tposted[0] = get_cycles();
	// main loop for posting 
	while (totscnt < (user_param->iters * user_param->num_of_qps)  || totccnt < (user_param->iters * user_param->num_of_qps) ) {
	    int iter;

	    // only poll cq if sent all
        if (totscnt < (user_param->iters * user_param->num_of_qps)) {

            // floating point exception, but if inner loop qps, performance down (lock increased)

    #ifdef _OPENMP_ON
    #pragma omp parallel for reduction(+: errs) shared(wr, rem_addr, my_addr) private(iter, index)
    #endif
            for (index = 0; index < user_param->num_of_qps; index++) {
                for (iter = 0; iter < iter_depth[index]; iter++) {
                    struct ibv_send_wr _wr;
                    struct ibv_send_wr *bad_wr;
                    int dpth = iter_depth[index];
                    _wr = wr[index];
                    _wr.wr_id = index + 1;

                    _wr.send_flags |= IBV_SEND_SIGNALED;

                    if (ibv_post_send(ctx->qp[index], &_wr, &bad_wr)) {
                        fprintf(stderr,
                                "%d)Couldn't post send: (mix:%d, depth:%d/%d, qp:%d), "
                                        "wr num_seg %d, addr 0x%lx, length %d, lkey 0x%x, flag 0x%x \n",
                                main_iter, iter, dpth, iter_depth[index], index,
                                _wr.num_sge, _wr.sg_list[0].addr,
                                _wr.sg_list[0].length, _wr.sg_list[0].lkey,
                                _wr.send_flags);
                        errs++;
                    } else if (_wr.wr_id != index + 1) {
                        fprintf(stderr, "2)wr_id id wrong %ld != %d\n",
                                _wr.wr_id, index + 1);
                        errs++;
                    }
                }
            }
            if (errs > 0) {
                return 1;
            }

            debug_printf("%d-1) totscnt %d totccnt %d\n", main_iter, totscnt, totccnt);
            // set qp full
            for (index = 0; index < user_param->num_of_qps; index++) {
                debug_printf("depth[%d] = %d \n", index, iter_depth[index]);
                totscnt += iter_depth[index];
                iter_depth[index] = 0;
            }
        }
		// clear cq
        do {
            ne = ibv_poll_cq(ctx->send_cq, DEF_WC_SIZE, wc);
            if (ne > 0) {
                for (i = 0; i < ne; i++) {

                    if (wc[i].status != IBV_WC_SUCCESS)
                        NOTIFY_COMP_ERROR_SEND(wc[i], totscnt, totccnt);

                    iter_depth[(int)wc[i].wr_id - 1]++; // set empty
                    totccnt++;
                }
            } else if (ne < 0) {
                break;
            }
        } while (totccnt < totscnt);
        if (ne < 0) {
            fprintf(stderr, "poll CQ failed %d\n", ne);
            return 1;
        }

        debug_printf("%d-2) totscnt %d totccnt %d\n", main_iter, totscnt, totccnt);
#ifdef DEBUG
        for (index = 0; index < user_param->num_of_qps; index++) {
            debug_printf("depth[%d] = %d \n", index, iter_depth[index]);
        }
#endif
        main_iter++;
	}
	tcompleted[user_param->iters*user_param->num_of_qps - 1] = get_cycles();

	free(wr);
	free(sge_list);
	free(my_addr);
	free(rem_addr);
	free(wc);
	free(iter_depth);
	return 0;
}

/****************************************************************************** 
 *
 ******************************************************************************/
int main(int argc, char *argv[]) {

	int                         i = 0;
	struct ibv_device		    *ib_dev = NULL;
	struct pingpong_context     ctx;
	struct pingpong_dest        *my_dest,*rem_dest;
	struct perftest_parameters  user_param;
	struct perftest_comm		user_comm;

	/* init default values to user's parameters */
	memset(&user_param,0,sizeof(struct perftest_parameters));
	memset(&user_comm,0,sizeof(struct perftest_comm));
	memset(&ctx,0,sizeof(struct pingpong_context));
	
	user_param.verb    = WRITE;
	user_param.tst     = BW;
	user_param.version = VERSION;

	// Configure the parameters values according to user arguments or defalut values.
	if (parser(&user_param,argv,argc)) {
		fprintf(stderr," Parser function exited with Error\n");
		return 1;
	}

	// Finding the IB device selected (or defalut if no selected).
	ib_dev = ctx_find_dev(user_param.ib_devname);
	if (!ib_dev) {
		fprintf(stderr," Unable to find the Infiniband/RoCE deivce\n");
		return 1;
	}

	// Getting the relevant context from the device
	ctx.context = ibv_open_device(ib_dev);
	if (!ctx.context) {
		fprintf(stderr, " Couldn't get context for the device\n");
		return 1;
	}

	// See if MTU and link type are valid and supported.
	if (check_link_and_mtu(ctx.context,&user_param)) {
		fprintf(stderr, " Couldn't get context for the device\n");
		return FAILURE;
	}

	// Print basic test information.
//	ctx_print_test_info(&user_param);

	ALLOCATE(my_dest , struct pingpong_dest , user_param.num_of_qps);
	memset(my_dest, 0, sizeof(struct pingpong_dest)*user_param.num_of_qps);
	ALLOCATE(rem_dest , struct pingpong_dest , user_param.num_of_qps);
	memset(rem_dest, 0, sizeof(struct pingpong_dest)*user_param.num_of_qps);

	// copy the rellevant user parameters to the comm struct + creating rdma_cm resources.
	if (create_comm_struct(&user_comm,&user_param)) { 
		fprintf(stderr," Unable to create RDMA_CM resources\n");
		return 1;
	}

	// Create (if nessacery) the rdma_cm ids and channel.
	if (user_param.work_rdma_cm == ON) {

	    if (create_rdma_resources(&ctx,&user_param)) {
			fprintf(stderr," Unable to create the rdma_resources\n");
			return FAILURE;
	    }
		
  	    if (user_param.machine == CLIENT) {

			if (rdma_client_connect(&ctx,&user_param)) {
				fprintf(stderr,"Unable to perform rdma_client function\n");
				return FAILURE;
			}
		
		} else {

			if (rdma_server_connect(&ctx,&user_param)) {
				fprintf(stderr,"Unable to perform rdma_client function\n");
				return FAILURE;
			}
		}
					
	} else {
    
	    // create all the basic IB resources (data buffer, PD, MR, CQ and events channel)
	    if (ctx_init(&ctx,&user_param)) {
			fprintf(stderr, " Couldn't create IB resources\n");
			return FAILURE;
	    }
	}

	// Set up the Connection.
	if (set_up_connection(&ctx,&user_param,my_dest)) {
		fprintf(stderr," Unable to set up socket connection\n");
		return FAILURE;
	}

	// Print this machine QP information
	for (i=0; i < user_param.num_of_qps; i++) 
		ctx_print_pingpong_data(&my_dest[i],&user_comm);

	// Init the connection and print the local data.
	if (establish_connection(&user_comm)) {
		fprintf(stderr," Unable to init the socket connection\n");
		return FAILURE;
	}

	// shaking hands and gather the other side info.
	for (i=0; i < user_param.num_of_qps; i++) {

			if (ctx_hand_shake(&user_comm,&my_dest[i],&rem_dest[i])) {
				fprintf(stderr," Failed to exchange date between server and clients\n");
				return 1;   
			}

			// Print remote machine QP information
			user_comm.rdma_params->side = REMOTE;
			ctx_print_pingpong_data(&rem_dest[i],&user_comm);

			if (user_param.work_rdma_cm == OFF) {

				if (pp_connect_ctx(&ctx,my_dest[i].psn,&rem_dest[i],&user_param,i)) {
					fprintf(stderr," Unable to Connect the HCA's through the link\n");
					return FAILURE;
				}
			}

			// An additional handshake is required after moving qp to RTR.
			if (ctx_hand_shake(&user_comm,&my_dest[i],&rem_dest[i])) {
				fprintf(stderr," Failed to exchange date between server and clients\n");
				return FAILURE; 
			}
	}	

	printf(RESULT_LINE);
	printf(RESULT_FMT);

	// For half duplex tests, server just waits for client to exit 
	if (user_param.machine == SERVER && !user_param.duplex) {
		
		if (ctx_close_connection(&user_comm,&my_dest[0],&rem_dest[0])) {
			fprintf(stderr,"Failed to close connection between server and client\n");
			return 1;
		}
		printf(RESULT_LINE);
		return destroy_ctx(&ctx,&user_param);
	}

	ALLOCATE(tposted,cycles_t,user_param.iters*user_param.num_of_qps);
	ALLOCATE(tcompleted,cycles_t,user_param.iters*user_param.num_of_qps);

#ifdef _OPENMP_ON
            int nth;
            for(nth = 1; nth < 240; nth*=2){
                omp_set_num_threads(nth);
#endif

    if (user_param.all == ON) {

        for (i = 1; i < 12; ++i) {
            user_param.size = 1 << i;
            if (run_iter(&ctx, &user_param, rem_dest))
                return 17;
            print_report(&user_param);
        }

    } else {

        if (run_iter(&ctx, &user_param, rem_dest))
            return 18;
        print_report(&user_param);
    }

#ifdef _OPENMP_ON
            }

            omp_set_num_threads(240);
            if (user_param.all == ON) {

                for (i = 1; i < 12 ; ++i) {
                    user_param.size = 1 << i;
                    if(run_iter(&ctx,&user_param,rem_dest))
                        return 17;
                    print_report(&user_param);
                }

            } else {

                if(run_iter(&ctx,&user_param,rem_dest))
                    return 18;
                print_report(&user_param);
            }
#endif

	free(tposted);
	free(tcompleted);

	// Closing connection.
	if (ctx_close_connection(&user_comm,&my_dest[0],&rem_dest[0])) {
	 	fprintf(stderr,"Failed to close connection between server and client\n");
		return 1;
	}

	free(my_dest);
	free(rem_dest);
	printf(RESULT_LINE);
	return destroy_ctx(&ctx,&user_param);
}

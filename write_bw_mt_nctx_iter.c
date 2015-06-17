/*
 * Copyright (c) 2005 Topspin Communications.  All rights reserved.
 * Copyright (c) 2005 Mellanox Technologies Ltd.  All rights reserved.
 * Copyright (c) 2009 HNR Consulting.  All rights reserved.
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
 * This experiment proves that multi-threading can improve BW for small messages.
 *
 * This test parallels ibv_post_send in ib context level.
 * Every thread has separate ibv_context/qp/cq, and then issues ibv_post_send in parallel.
 */

#if HAVE_CONFIG_H
#  include <config.h>
#endif /* HAVE_CONFIG_H */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <malloc.h>

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
#define NCTX 64

cycles_t    *tposted;
cycles_t    *tcompleted;

// use separate cq for each qp
static struct ibv_cq **nsend_cq;
int g_totscnt; // the real sent time

/****************************************************************************** 
 *
 ******************************************************************************/
static int pp_connect_ctx(struct pingpong_context *ctx, int my_psn,
        struct pingpong_dest *dest, struct perftest_parameters *user_parm,
        int qpindex) {
    struct ibv_qp_attr attr;
    memset(&attr, 0, sizeof attr);

    attr.qp_state = IBV_QPS_RTR;
    attr.path_mtu = user_parm->curr_mtu;
    attr.dest_qp_num = dest->qpn;
    attr.rq_psn = dest->psn;
    attr.ah_attr.dlid = dest->lid;
    if (user_parm->connection_type == RC) {
        attr.max_dest_rd_atomic = 1;
        attr.min_rnr_timer = 12;
    }
    if (user_parm->gid_index < 0) {
        attr.ah_attr.is_global = 0;
        attr.ah_attr.sl = user_parm->sl;
    } else {
        attr.ah_attr.is_global = 1;
        attr.ah_attr.grh.dgid = dest->gid;
        attr.ah_attr.grh.sgid_index = user_parm->gid_index;
        attr.ah_attr.grh.hop_limit = 1;
        attr.ah_attr.sl = 0;
    }
    attr.ah_attr.src_path_bits = 0;
    attr.ah_attr.port_num = user_parm->ib_port;
    if (user_parm->connection_type == RC) {
        if (ibv_modify_qp(ctx->qp[qpindex], &attr,
                IBV_QP_STATE | IBV_QP_AV | IBV_QP_PATH_MTU | IBV_QP_DEST_QPN
                        | IBV_QP_RQ_PSN | IBV_QP_MIN_RNR_TIMER
                        | IBV_QP_MAX_DEST_RD_ATOMIC)) {
            fprintf(stderr, "Failed to modify RC QP to RTR\n");
            return 1;
        }
        attr.timeout = user_parm->qp_timeout;
        attr.retry_cnt = 7;
        attr.rnr_retry = 7;
    } else {
        if (ibv_modify_qp(ctx->qp[qpindex], &attr,
                IBV_QP_STATE | IBV_QP_AV | IBV_QP_PATH_MTU | IBV_QP_DEST_QPN
                        | IBV_QP_RQ_PSN)) {
            fprintf(stderr, "Failed to modify UC QP to RTR\n");
            return 1;
        }

    }
    attr.qp_state = IBV_QPS_RTS;
    attr.sq_psn = my_psn;
    attr.max_rd_atomic = 1;
    if (user_parm->connection_type == 0) {
        attr.max_rd_atomic = 1;
        if (ibv_modify_qp(ctx->qp[qpindex], &attr,
                IBV_QP_STATE | IBV_QP_SQ_PSN | IBV_QP_TIMEOUT | IBV_QP_RETRY_CNT
                        | IBV_QP_RNR_RETRY | IBV_QP_MAX_QP_RD_ATOMIC)) {
            fprintf(stderr, "Failed to modify RC QP to RTS\n");
            return 1;
        }
    } else {
        if (ibv_modify_qp(ctx->qp[qpindex], &attr,
                IBV_QP_STATE | IBV_QP_SQ_PSN)) {
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
    unsigned long tsize; /* Transferred size, in megabytes */
    int i, j;
    int opt_posted = 0, opt_completed = 0;
    cycles_t opt_delta;
    cycles_t t;
    int iters = user_param->iters;
    int nthreads;
    int exceedcnt;

#pragma omp parallel
    nthreads = omp_get_num_threads();

    cycles_to_units = get_cpu_mhz(user_param->cpu_freq_f) * 1000000;

    tsize = user_param->duplex ? 2 : 1;
    tsize = tsize * user_param->size;

    int time = tcompleted[1] - tposted[1];
    int niter = g_totscnt - user_param->tx_depth * user_param->num_of_qps * NCTX;
    double bw = tsize * niter * cycles_to_units / time / 0x100000;

    printf("%d %lu %d %.2f\n", nthreads, (unsigned long) user_param->size, iters, bw);
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
inline void increase_loc_addr_ncnt(struct ibv_sge *sg, int size, int rcnt,
        uint64_t prim_addr, int server_is_ud) {

    sg->addr = prim_addr + INC(size) * rcnt % CYCLE_BUFFER;
}


#undef DEF_WC_SIZE
#define DEF_WC_SIZE         (300)

/******************************************************************************
 *
 ******************************************************************************/
int run_iter(struct pingpong_context ctx[NCTX],
             struct perftest_parameters *user_param,
             struct pingpong_dest *rem_dest)    {

    int                totscnt = 0;
    int                totccnt = 0;
    int                i       = 0, j = 0;
    int                index,ne;
    int                warmindex;
    struct ibv_send_wr *bad_wr;
    struct ibv_wc      wc[NCTX * DEF_WC_SIZE];
    struct ibv_sge     *sge_list = NULL;
    struct ibv_send_wr *wr       = NULL;
    uint64_t           *my_addr  = NULL;
    uint64_t           *rem_addr = NULL;
    int *iter_depth; /* available depth of each iteration */
    int errs = 0;

    ALLOCATE(wr, struct ibv_send_wr, NCTX*user_param->num_of_qps);
    ALLOCATE(sge_list, struct ibv_sge, NCTX*user_param->num_of_qps);
    ALLOCATE(my_addr, uint64_t, NCTX*user_param->num_of_qps);
    ALLOCATE(rem_addr, uint64_t, NCTX*user_param->num_of_qps);
    ALLOCATE(iter_depth, int, user_param->num_of_qps*NCTX);
    memset(wc, 0, sizeof(struct ibv_wc) * NCTX * DEF_WC_SIZE);


    // Each CTX/QP has its own wr and sge , that holds the qp addresses and attr.
    // We write in cycles on the buffer to exploid the "Nahalem" system.
    for (i = 0; i < NCTX; i++) {
        for (j = 0; j < user_param->num_of_qps; j++) {
            index = i * user_param->num_of_qps + j;

            sge_list[index].addr = (uintptr_t) ctx[i].buf
                    + (j * BUFF_SIZE(ctx[i].size));
            sge_list[index].length = user_param->size;
            sge_list[index].lkey = ctx[i].mr->lkey;

            wr[index].sg_list = &sge_list[index];
            wr[index].num_sge = 1;
            wr[index].opcode = IBV_WR_RDMA_WRITE;
            wr[index].next = NULL;
            wr[index].wr.rdma.remote_addr = rem_dest[index].vaddr;
            wr[index].wr.rdma.rkey = rem_dest[index].rkey;
            wr[index].wr_id = index + 1;
            wr[index].send_flags = IBV_SEND_SIGNALED;

            if (user_param->size <= user_param->inline_size)
                wr[index].send_flags |= IBV_SEND_INLINE;

            debug_printf("prepare: ctx:%d qp:%d, wr_id %d, addr 0x%llx, len %d, lkey 0x%x, rdma_addr 0x%llx, rkey 0x%x\n"
                    , i, j, wr[index].wr_id,
                    sge_list[index].addr, sge_list[index].length, sge_list[index].lkey,
                    wr[index].wr.rdma.remote_addr, wr[index].wr.rdma.rkey);

            iter_depth[index] = user_param->tx_depth;
        }
    }

    tposted[0] = get_cycles();
#ifdef _OPENMP_ON
#pragma omp parallel for reduction(+: errs) shared(ctx, wr, rem_addr, my_addr, iter_depth, user_param)  \
    private(warmindex, index, j)
#endif
    for (i = 0; i < NCTX; i++) {
        for (j = 0; j < user_param->num_of_qps; j++) {
            index = i * user_param->num_of_qps + j;

            for (warmindex = 0; warmindex < user_param->tx_depth; warmindex++) {
                /* private per thread */
                struct ibv_send_wr _wr = wr[index];
                struct ibv_send_wr *bad_wr;
                _wr.send_flags |= IBV_SEND_SIGNALED;
                _wr.wr_id = index + 1;

                if (ibv_post_send(ctx[i].qp[j], &_wr, &bad_wr)) {
                    fprintf(stderr, "1)[%d/%d]Couldn't post send: (depth:%d, ctx:%d qp:%d), "
                            "wr num_seg %d, addr 0x%lx, length %d, lkey 0x%x, flag 0x%x \n",
                            omp_get_thread_num(), omp_get_num_threads(),
                            warmindex, i, j, _wr.num_sge, _wr.sg_list[0].addr,
                            _wr.sg_list[0].length, _wr.sg_list[0].lkey, _wr.send_flags);
                    errs++;
                    break;
                }
//                debug_printf("1)[%d/%d]post send: (depth:%d, ctx:%d qp:%d), "
//                        "wr num_seg %d, addr 0x%lx, length %d, lkey 0x%x, flag 0x%x, rdma_addr 0x%llx, rkey 0x%x wr_id %d\n",
//                        omp_get_thread_num(), omp_get_num_threads(),
//                        warmindex, i, j, _wr.num_sge, _wr.sg_list[0].addr,
//                        _wr.sg_list[0].length, _wr.sg_list[0].lkey, _wr.send_flags,
//                        _wr.wr.rdma.remote_addr, _wr.wr.rdma.rkey, _wr.wr_id);
            }
            iter_depth[index] -= user_param->tx_depth;

            debug_printf("[%d/%d] warm: ctx:%d %p qp:%d %p, depth %d/%d\n",
                    omp_get_thread_num(), omp_get_num_threads(),
                    i, &ctx[i], j, &ctx[i].qp[j], iter_depth[index] , user_param->tx_depth);
        }
    }
    if (errs > 0) {
        fprintf(stderr, "warm posting send failed, errs > 0\n");
        return 1;
    }

    totscnt = user_param->tx_depth * user_param->num_of_qps * NCTX;

    int xx = -1;
    // clear cq
    while (totccnt < totscnt && errs == 0) {
        for (i = 0; i < NCTX && errs == 0; i++) {
            int wc_start = i * DEF_WC_SIZE;
            ne = ibv_poll_cq(ctx[i].send_cq, DEF_WC_SIZE, &wc[wc_start]);
            if (ne > 0) {
                for ( j= 0; j < ne; j++) {
                    xx++;
                    if (wc[wc_start + j].status != IBV_WC_SUCCESS){
                        fprintf(stderr, "CTX[%d] cqe[%d] wc_id %ld status wrong, %d\n", i, xx,
                                wc[wc_start + j].wr_id,
                                wc[wc_start + j].status);
                        errs++;
                        break;
                    }
//                    debug_printf("warm: ctx:%d, polled cqe[%d] wc[%d] wr_id %d\n", i, xx,
//                            wc_start + i, wc[wc_start + i].wr_id);
                }
                totccnt +=ne;
                iter_depth[i] += ne;
            }
            // ne < 0, wrong
            else if (ne < 0) {
                errs++;
                fprintf(stderr, "CTX[%d] poll CQ failed, error %d\n", i, ne);
                break;
            }
        }
    }
    if (errs > 0) {
        fprintf(stderr, "poll CQ failed, errs %d > 0\n", errs);
        return 1;
    }

    tcompleted[0] = get_cycles();

    debug_printf("1-2) totscnt %d totccnt %d\n", totscnt, totccnt);
#pragma omp parallel for shared(iter_depth, user_param) private(index)
    for (i = 0; i < NCTX; i++) {
        for (j = 0; j < user_param->num_of_qps; j++) {
            index = i * user_param->num_of_qps + j;
            iter_depth[index] = user_param->tx_depth;
        }
    }

    int main_iter = 2;
    tposted[1] = get_cycles();
    // main loop for posting
    while (totscnt < (user_param->iters * user_param->num_of_qps * NCTX)  ||
            totccnt < (user_param->iters * user_param->num_of_qps * NCTX) ||
            totccnt < totscnt /*totscnt exceed sometimes*/) {
        int iter;

        // only poll cq if sent all
        if (totscnt < (user_param->iters * user_param->num_of_qps * NCTX)) {

#ifdef _OPENMP_ON
#pragma omp parallel for reduction(+: errs, totscnt, totccnt) \
            shared(wr, rem_addr, my_addr, nsend_cq, iter_depth, user_param, wc) private(iter, index, j)
#endif
            for (i = 0; i < NCTX; i++) {
                for (j = 0; j < user_param->num_of_qps; j++) {
                    index = i * user_param->num_of_qps + j;

                    for (iter = 0; iter < iter_depth[index]; iter++) {
                        struct ibv_send_wr _wr = wr[index];
                        struct ibv_send_wr *bad_wr;
                        int dpth = iter_depth[index];
                        _wr.wr_id = iter + 1;
                        _wr.send_flags |= IBV_SEND_SIGNALED;

//                        if ((iter+1) % user_param->cq_mod == 0 && user_param->cq_mod > 1){
//                            _wr.send_flags &= IBV_SEND_SIGNALED;
//                            printf( "%d)[%d/%d] ctx:%d qp:%d post SIGNALED wr_id %d\n",
//                                   main_iter, omp_get_thread_num(), omp_get_num_threads(), i, j,
//                                   _wr.wr_id);
//                        }

                        if (ibv_post_send(ctx[i].qp[j], &_wr, &bad_wr)) {
                            fprintf(stderr,
                                    "%d)[%d/%d]Couldn't post send: (mix:%d, depth:%d/%d, ctx:%d qp:%d), "
                                            "wr num_seg %d, addr 0x%lx, length %d, lkey 0x%x, flag 0x%x \n",
                                    main_iter, omp_get_thread_num(),
                                    omp_get_num_threads(), iter, dpth,
                                    iter_depth[index], i, j, _wr.num_sge,
                                    _wr.sg_list[0].addr, _wr.sg_list[0].length,
                                    _wr.sg_list[0].lkey, _wr.send_flags);
                            errs++;
                            break;
                        }
                    }

                    // set qp full
                    totscnt += iter_depth[index];
                    iter_depth[index] = 0;

                    debug_printf(
                            "%d-1-1)[%d/%d] ctx:%d qp:%d, depth %d/%d\n", main_iter, omp_get_thread_num(), omp_get_num_threads(), i, j, iter_depth[index], user_param->tx_depth);
                }

                do {
                    int wc_start = i * DEF_WC_SIZE;
                    ne = ibv_poll_cq(ctx[i].send_cq, DEF_WC_SIZE,
                            &wc[wc_start]);
                    if (ne > 0) {
                        for (j = 0; j < ne; j++) {
                            if (wc[wc_start + j].status != IBV_WC_SUCCESS) {
                                fprintf(stderr, "CTX[%d] cq status wrong, %d\n",
                                        j, wc[wc_start + j].status);
                                errs++;
                                break;
                            }
                            //                            printf( "CTX[%d] cq polled wr_id %d\n",i, wc[wc_start + i].wr_id);
                        }
                        iter_depth[i] += ne;
                        totccnt += ne;
                    }
                    // ne < 0, wrong
                    else if (ne < 0) {
                        errs++;
                        fprintf(stderr, "CTX[%d] poll CQ failed, error %d\n", i,
                                ne);
                        break;
                    }
                } while (ne > 0);

                debug_printf(
                        "%d-1-2)[%d/%d] ctx:%d, depth %d/%d\n", main_iter, omp_get_thread_num(), omp_get_num_threads(), i, iter_depth[index], user_param->tx_depth);
            }
            if (errs > 0) {
                fprintf(stderr, "posting send failed, errs > 0\n");
                return 1;
            }
            debug_printf(
                    "%d-1) totscnt %d totccnt %d\n", main_iter, totscnt, totccnt);
        }

        // clear cq as many as we can
        // todo: it clears all, need to change, but may ctx[0] can only free 1 if do not
        do {
#ifdef _OPENMP_ON
//#pragma omp parallel for reduction(+: errs, totccnt) \
//        shared(wr, rem_addr, my_addr, nsend_cq, iter_depth, user_param, wc) private(iter, index)
#endif
            for (i = 0; i < NCTX; i++) {
                int wc_start = i * DEF_WC_SIZE;
                if (iter_depth[i] == user_param->tx_depth) continue;
                do {
                    ne = ibv_poll_cq(ctx[i].send_cq, DEF_WC_SIZE,
                            &wc[wc_start]);
                    if (ne > 0) {
                        for (j = 0; j < ne; j++) {
                            if (wc[wc_start + j].status != IBV_WC_SUCCESS) {
                                fprintf(stderr, "CTX[%d] cq status wrong, %d\n",
                                        j, wc[wc_start + j].status);
                                errs++;
                                break;
                            }
//                            printf( "CTX[%d] cq polled wr_id %d\n",i, wc[wc_start + i].wr_id);
                        }
                        iter_depth[i] += ne;
                        totccnt += ne;
                    }
                    // ne < 0, wrong
                    else if (ne < 0) {
                        errs++;
                        fprintf(stderr, "CTX[%d] poll CQ failed, error %d\n", i,
                                ne);
                        break;
                    }
                } while (ne > 0);
            }
        } while (errs == 0 && totscnt > totccnt);
        if (errs > 0) {
            fprintf(stderr, "poll CQ failed, errs > 0\n");
            return 1;
        }

        debug_printf("%d-2) totscnt %d totccnt %d\n", main_iter, totscnt, totccnt);
#ifdef DEBUG
        for (i = 0; i < NCTX; i++) {
            for (j = 0; j < user_param->num_of_qps; j++) {
                index = i * user_param->num_of_qps + j;
                debug_printf("ctx[%d].qp[%d] depth = %d \n", i, j, iter_depth[index]);
            }
        }
#endif
        main_iter++;
    }

    debug_printf("x) totscnt %d totccnt %d\n", totscnt, totccnt);
    tcompleted[1] = get_cycles();

    g_totscnt = totscnt;

    debug_printf("x) free wr\n");
    free(wr);
    debug_printf("x) free sge\n");
    free(sge_list);
    debug_printf("x) free my_addr\n");
    free(my_addr);
    debug_printf("x) free rem_addr\n");
    free(rem_addr);
    debug_printf("x) free iter_depth\n");
    free(iter_depth);
    debug_printf("x) free all\n");

    return 0;
}

/******************************************************************************
 *
 ******************************************************************************/
struct ibv_qp* ncq_ctx_qp_create(struct pingpong_context *ctx,
        struct ibv_cq* scq, struct perftest_parameters *user_param) {

    struct ibv_qp_init_attr attr;
    struct ibv_qp* qp = NULL;

    memset(&attr, 0, sizeof(struct ibv_qp_init_attr));
    attr.send_cq = scq;
    attr.recv_cq = (user_param->verb == SEND) ? ctx->recv_cq : scq;
    attr.cap.max_send_wr = user_param->tx_depth;
    attr.cap.max_recv_wr = user_param->rx_depth;
    attr.cap.max_send_sge = MAX_SEND_SGE;
    attr.cap.max_recv_sge = MAX_RECV_SGE;
    attr.cap.max_inline_data = user_param->inline_size;

    switch (user_param->connection_type) {

    case RC:
        attr.qp_type = IBV_QPT_RC;
        break;
    case UC:
        attr.qp_type = IBV_QPT_UC;
        break;
    case UD:
        attr.qp_type = IBV_QPT_UD;
        break;
    default:
        fprintf(stderr, "Unknown connection type \n");
        return NULL;
    }

    if (user_param->work_rdma_cm) {

        if (rdma_create_qp(ctx->cm_id, ctx->pd, &attr)) {
            fprintf(stderr, " Couldn't create rdma QP\n");
            perror("rdma_create_qp");
            return NULL;
        }

        qp = ctx->cm_id->qp;

    } else {

        qp = ibv_create_qp(ctx->pd, &attr);
    }
    return qp;
}

/****************************************************************************** 
 *
 ******************************************************************************/
int ncq_ctx_init(struct pingpong_context *ctx,
        struct perftest_parameters *user_param) {

    int i, flags;
    uint64_t buff_size;

    ALLOCATE(ctx->qp, struct ibv_qp*, user_param->num_of_qps);
    ALLOCATE(ctx->scnt, int, user_param->num_of_qps);
    ALLOCATE(ctx->ccnt, int, user_param->num_of_qps);
    ALLOCATE(nsend_cq, struct ibv_cq*, user_param->num_of_qps);

    memset(ctx->scnt, 0, user_param->num_of_qps * sizeof(int));
    memset(ctx->ccnt, 0, user_param->num_of_qps * sizeof(int));

    flags = IBV_ACCESS_REMOTE_WRITE | IBV_ACCESS_LOCAL_WRITE;
    // ctx->size = SIZE(user_param->connection_type,user_param->size,!(int)user_param->machine);
    ctx->size = user_param->size;
    buff_size = BUFF_SIZE(ctx->size) * 2 * user_param->num_of_qps;

    // Allocating the buffer in BUFF_SIZE size to support max performance.
    ctx->buf = memalign(sysconf(_SC_PAGESIZE), buff_size);
    if (!ctx->buf) {
        fprintf(stderr, "Couldn't allocate work buf.\n");
        return FAILURE;
    }
    memset(ctx->buf, 0, buff_size);

    // Allocating an event channel if requested.
    if (user_param->use_event) {
        ctx->channel = ibv_create_comp_channel(ctx->context);
        if (!ctx->channel) {
            fprintf(stderr, "Couldn't create completion channel\n");
            return FAILURE;
        }
    }

    // Allocating the Protection domain.
    ctx->pd = ibv_alloc_pd(ctx->context);
    if (!ctx->pd) {
        fprintf(stderr, "Couldn't allocate PD\n");
        return FAILURE;
    }

    if (user_param->verb == READ) {
        flags |= IBV_ACCESS_REMOTE_READ;

    } else if (user_param->verb == ATOMIC) {
        flags |= IBV_ACCESS_REMOTE_ATOMIC;
    }

    // Alocating Memory region and assiging our buffer to it.
    ctx->mr = ibv_reg_mr(ctx->pd, ctx->buf, buff_size, flags);
    if (!ctx->mr) {
        fprintf(stderr, "Couldn't allocate MR\n");
        return FAILURE;
    }

//    ctx->send_cq = ibv_create_cq(ctx->context,user_param->tx_depth*user_param->num_of_qps,NULL,ctx->channel,0);
//    if (!ctx->send_cq) {
//        fprintf(stderr, "Couldn't create CQ\n");
//        return FAILURE;
//    }

    for (i = 0; i < user_param->num_of_qps; i++) {
        nsend_cq[i] = ibv_create_cq(ctx->context, user_param->tx_depth, NULL,
                ctx->channel, 0);
        if (!nsend_cq[i]) {
            fprintf(stderr, "Couldn't create CQ[%d]\n", i);
            return FAILURE;
        }

        debug_printf("CQ[%d] %p created\n", i, nsend_cq[i]);
    }

    if (user_param->verb == SEND) {
        ctx->recv_cq = ibv_create_cq(ctx->context,
                user_param->rx_depth * user_param->num_of_qps, NULL,
                ctx->channel, 0);
        if (!ctx->recv_cq) {
            fprintf(stderr, "Couldn't create a recevier CQ\n");
            return FAILURE;
        }
    }

    for (i = 0; i < user_param->num_of_qps; i++) {

        ctx->qp[i] = ncq_ctx_qp_create(ctx, nsend_cq[i], user_param);
        if (ctx->qp[i] == NULL) {
            fprintf(stderr, " Unable to create QP.\n");
            return FAILURE;
        } debug_printf("QP[%d] %p created\n", i, ctx->qp[i]);

        if (user_param->work_rdma_cm == OFF) {

            if (ctx_modify_qp_to_init(ctx->qp[i], user_param)) {
                fprintf(stderr, "Failed to modify QP to INIT\n");
                return FAILURE;
            }
        }
    }

    return SUCCESS;
}

int ncq_destroy_ctx(struct pingpong_context *ctx,
        struct perftest_parameters *user_parm) {

    int i;
    int test_result = 0;

    for (i = 0; i < user_parm->num_of_qps; i++) {
        if (ibv_destroy_qp(ctx->qp[i])) {
            fprintf(stderr, "failed to destroy QP\n");
            test_result = 1;
        } debug_printf("QP[%d] %p destroyed\n", i, ctx->qp[i]);
        if (ibv_destroy_cq(nsend_cq[i])) {
            fprintf(stderr, "failed to destroy CQ %d\n", i);
            test_result = 1;
        } debug_printf("CQ[%d] %p created\n", i, nsend_cq[i]);
    }

//    if (ibv_destroy_cq(ctx->send_cq)) {
//        fprintf(stderr, "failed to destroy CQ\n");
//        test_result = 1;
//    }

    if (user_parm->verb == SEND) {
        if (ibv_destroy_cq(ctx->recv_cq)) {
            fprintf(stderr, "failed to destroy CQ\n");
            test_result = 1;
        }
    }

    if (ibv_dereg_mr(ctx->mr)) {
        fprintf(stderr, "failed to deregister MR\n");
        test_result = 1;
    }

    if (ibv_dealloc_pd(ctx->pd)) {
        fprintf(stderr, "failed to deallocate PD\n");
        test_result = 1;
    }

    if (ctx->channel) {
        if (ibv_destroy_comp_channel(ctx->channel)) {
            test_result = 1;
        }
    }

    if (user_parm->work_rdma_cm == OFF) {

        if (ibv_close_device(ctx->context)) {
            fprintf(stderr, "failed to close device context\n");
            test_result = 1;
        }

    } else {

        rdma_destroy_id(ctx->cm_id);
        rdma_destroy_event_channel(ctx->cm_channel);
    }

    free(ctx->buf);
    free(ctx->qp);
    free(ctx->scnt);
    free(ctx->ccnt);

    if (nsend_cq) {
        free(nsend_cq);
    }

    return test_result;
}
/******************************************************************************
 *
 ******************************************************************************/

int main(int argc, char *argv[]) {

    int i = 0, j = 0;
    struct ibv_device *ib_dev = NULL;
    struct pingpong_context ctx[NCTX];
    struct pingpong_dest my_dest[NCTX], rem_dest[NCTX];
    struct perftest_parameters user_param;
    struct perftest_comm user_comm;

    /* init default values to user's parameters */
    memset(&user_param, 0, sizeof(struct perftest_parameters));
    memset(&user_comm, 0, sizeof(struct perftest_comm));
    memset(ctx, 0, sizeof(struct pingpong_context) * NCTX);

    user_param.verb = WRITE;
    user_param.tst = BW;
    user_param.version = VERSION;

    // Configure the parameters values according to user arguments or defalut values.
    if (parser(&user_param, argv, argc)) {
        fprintf(stderr, " Parser function exited with Error\n");
        return 1;
    }

    user_param.num_of_qps = 1;
//    user_param.tx_depth = 4;
//    user_param.cq_mod = 1;

    // Finding the IB device selected (or defalut if no selected).
    ib_dev = ctx_find_dev(user_param.ib_devname);
    if (!ib_dev) {
        fprintf(stderr, " Unable to find the Infiniband/RoCE deivce\n");
        return 1;
    }

    memset(my_dest, 0,
            sizeof(struct pingpong_dest) * user_param.num_of_qps * NCTX);
    memset(rem_dest, 0,
            sizeof(struct pingpong_dest) * user_param.num_of_qps * NCTX);

    for (i = 0; i < NCTX; i++) {
        // Getting the relevant context from the device
        ctx[i].context = ibv_open_device(ib_dev);
        if (!ctx[i].context) {
            fprintf(stderr, " Couldn't get context for the device\n");
            return 1;
        }

        // See if MTU and link type are valid and supported.
        if (check_link_and_mtu(ctx[i].context, &user_param)) {
            fprintf(stderr, " Couldn't get context for the device\n");
            return FAILURE;
        }

        // create all the basic IB resources (data buffer, PD, MR, CQ and events channel)
        if (ctx_init(&ctx[i], &user_param)) {
            fprintf(stderr, " Couldn't create IB resources\n");
            return FAILURE;
        }

        // Set up the Connection.
        if (set_up_connection(&ctx[i], &user_param, &my_dest[i])) {
            fprintf(stderr, " Unable to set up socket connection\n");
            return FAILURE;
        }

        debug_printf("init ctx[%d] context %p, pd %p, qp %p, cq %p/%p"
                ", mr %p (addr %p key 0x%x/0x%x), vaddr 0x%llx\n"
                , i, ctx[i].context, ctx[i].pd, ctx[i].qp[0], ctx[i].send_cq, ctx[i].qp[0]->send_cq
                ,ctx[i].mr, ctx[i].mr->addr, ctx[i].mr->lkey, ctx[i].mr->rkey, my_dest[i].vaddr
                );

        // Print basic test information.
//        ctx_print_test_info(&user_param);
    }

    // copy the rellevant user parameters to the comm struct + creating rdma_cm resources.
    if (create_comm_struct(&user_comm, &user_param)) {
        fprintf(stderr, " Unable to create RDMA_CM resources\n");
        return 1;
    }

    // Init the connection and print the local data.
    if (establish_connection(&user_comm)) {
        fprintf(stderr, " Unable to init the socket connection\n");
        return FAILURE;
    }

    for (i = 0; i < NCTX; i++) {
        for (j = 0; j < user_param.num_of_qps; j++) {
            int idx = i * user_param.num_of_qps + j;
            if (ctx_hand_shake(&user_comm, &my_dest[idx], &rem_dest[idx])) {
                fprintf(stderr,
                        " Failed to exchange date between server and clients\n");
                return 1;
            }

            // Print remote machine QP information
            user_comm.rdma_params->side = REMOTE;
            ctx_print_pingpong_data(&rem_dest[idx], &user_comm);

            if (user_param.work_rdma_cm == OFF) {
                if (pp_connect_ctx(&ctx[i], my_dest[idx].psn, &rem_dest[idx],
                        &user_param, j)) {
                    fprintf(stderr,
                            " Unable to Connect the HCA's through the link\n");
                    return FAILURE;
                }
            }

            // An additional handshake is required after moving qp to RTR.
            if (ctx_hand_shake(&user_comm, &my_dest[idx], &rem_dest[idx])) {
                fprintf(stderr,
                        " Failed to exchange date between server and clients\n");
                return FAILURE;
            }
        }
    }

    printf(RESULT_LINE);
    printf(RESULT_FMT);

    // For half duplex tests, server just waits for client to exit
    if (user_param.machine == SERVER && !user_param.duplex) {

        if (ctx_close_connection(&user_comm, &my_dest[0], &rem_dest[0])) {
            fprintf(stderr,
                    "Failed to close connection between server and client\n");
            return 1;
        }
        printf(RESULT_LINE);
        for(i = 0; i< NCTX; i++)
            destroy_ctx(&ctx[i], &user_param);
        return 0;
    }

    ALLOCATE(tposted, cycles_t, 2*NCTX);
    ALLOCATE(tcompleted, cycles_t, 2*NCTX);

    int ret;
#ifdef _OPENMP_ON
            int nth;
            for(nth = 1; nth < 128; nth*=2){
                omp_set_num_threads(nth);
#endif

    if (user_param.all == ON) {
        for (i = 1; i < 12; ++i) {
            user_param.size = 1 << i;
                ret = run_iter(ctx, &user_param, rem_dest);
            if (ret)
                return 17;
            print_report(&user_param);
        }
    } else {

        if (run_iter(ctx, &user_param, rem_dest))
            return 18;
        print_report(&user_param);
    }
#ifdef _OPENMP_ON
            }

            omp_set_num_threads(240);
            if (user_param.all == ON) {
                for (i = 1; i < 12; ++i) {
                    user_param.size = 1 << i;
                        ret = run_iter(ctx, &user_param, rem_dest);
                    if (ret)
                        return 17;
                    print_report(&user_param);
                }
            } else {

                if (run_iter(ctx, &user_param, rem_dest))
                    return 18;
                print_report(&user_param);
            }
#endif

    free(tposted);
    free(tcompleted);

    // Closing connection.
    if (ctx_close_connection(&user_comm, &my_dest[0], &rem_dest[0])) {
        fprintf(stderr,
                "Failed to close connection between server and client\n");
        return 1;
    }

    printf(RESULT_LINE);
    for(i = 0; i< NCTX; i++)
        destroy_ctx(&ctx[i], &user_param);
    return 0;
}

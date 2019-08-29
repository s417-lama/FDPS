#include "key.hpp"

#pragma once

namespace ParticleSimulator{

#ifdef PARTICLE_SIMULATOR_TASK_PARALLEL

    /* From ExaFMM's implemantation */

#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
#define N_TREE_CHILDREN 4
#else
#define N_TREE_CHILDREN 8
#endif

    struct BinaryTreeNode {
        S32              NBODY[N_TREE_CHILDREN];
        BinaryTreeNode * LEFT;
        BinaryTreeNode * RIGHT;
        BinaryTreeNode * BEGIN;
        BinaryTreeNode * END;
    };

    struct TreeNode {
        S32        IBODY;
        S32        NBODY;
        S32        NNODE;
        TreeNode * CHILD[N_TREE_CHILDREN];
    };

    inline S32 getNumBinNode(const S32 n) {
        if (n <= CUTOFF_SORT) return 1;
        else return 4 * ((n - 1) / CUTOFF_SORT) - 1;
    }

    template<class Tobj>
    void countBodies(const Tobj * val, const S32 begin, const S32 end,
                     const S32 level, BinaryTreeNode * binNode) {
        assert(getNumBinNode(end - begin) <= binNode->END - binNode->BEGIN + 1);
        if (end - begin <= CUTOFF_SORT) {
            for (S32 i = 0; i < N_TREE_CHILDREN; i++) binNode->NBODY[i] = 0;
            binNode->LEFT = binNode->RIGHT = NULL;
            for (S32 i = begin; i < end; i++) {
                S32 k = MortonKey::getCellID(level, val[i].getKey());
                assert(0 <= k && k < N_TREE_CHILDREN);
                binNode->NBODY[k]++;
            }
        } else {
            S32 mid = (begin + end) / 2;
            S32 numLeftNode  = getNumBinNode(mid - begin);
            S32 numRightNode = getNumBinNode(end - mid);
            assert(numLeftNode + numRightNode <= binNode->END - binNode->BEGIN);
            binNode->LEFT         = binNode->BEGIN;
            binNode->LEFT->BEGIN  = binNode->LEFT + 1;
            binNode->LEFT->END    = binNode->LEFT + numLeftNode;
            binNode->RIGHT        = binNode->LEFT->END;
            binNode->RIGHT->BEGIN = binNode->RIGHT + 1;
            binNode->RIGHT->END   = binNode->RIGHT + numRightNode;
            mk_task_group_w(2);
            create_taskc_w([&] { countBodies(val, begin, mid, level, binNode->LEFT ); }, 1);
            create_taskc_w([&] { countBodies(val, mid  , end, level, binNode->RIGHT); }, 1);
            wait_tasks;
            for (S32 i = 0; i < N_TREE_CHILDREN; i++) {
                binNode->NBODY[i] = binNode->LEFT->NBODY[i] + binNode->RIGHT->NBODY[i];
            }
        }
    }

    template<class Tobj>
    void moveBodies(const Tobj * val, Tobj * val_buf, const S32 begin, const S32 end,
                    const S32 level, BinaryTreeNode * binNode, const S32 * init_offset) {
        S32 offset[N_TREE_CHILDREN];
        for (S32 i = 0; i < N_TREE_CHILDREN; i++) {
            offset[i] = init_offset[i];
        }
        if (binNode->LEFT == NULL) {
            for (S32 i = begin; i < end; i++) {
                S32 k = MortonKey::getCellID(level, val[i].getKey());
                assert(0 <= k && k < N_TREE_CHILDREN);
                val_buf[offset[k]] = val[i];
                offset[k]++;
            }
        } else {
            int mid = (begin + end) / 2;
            mk_task_group_w(2);
            create_taskc_w([=] {
                moveBodies(val, val_buf, begin, mid, level, binNode->LEFT, offset);
            }, 1);
            for (S32 i = 0; i < N_TREE_CHILDREN; i++) offset[i] += binNode->LEFT->NBODY[i];
            create_taskc_w([=] {
                moveBodies(val, val_buf, mid, end, level, binNode->RIGHT, offset);
            }, 1);
            wait_tasks;
        }
    }

    inline void exclusiveScan(S32 * output, S32 * input, S32 offset) {
        for (S32 i = 0; i < N_TREE_CHILDREN; i++) {
            output[i] = offset;
            offset += input[i];
        }
    }

    inline S32 getMaxBinNode(const S32 n) {
        return (4 * n) / CUTOFF_SORT;
    }

    TreeNode * makeTreeNode(const bool nochild, const S32 begin, const S32 end) {
        TreeNode * node = new TreeNode();
        node->IBODY = begin;
        node->NBODY = end - begin;
        node->NNODE = 1;
        if (nochild) {
            for (S32 i = 0; i < N_TREE_CHILDREN; i++) node->CHILD[i] = NULL;
        }
        return node;
    }

    template<class Tobj>
    void sortBodiesImpl(Tobj * val, Tobj * val_buf, const S32 begin, const S32 end,
                        const S32 level, BinaryTreeNode * binNode, const S32 n_leaf_limit,
                        const bool direction, TreeNode ** node) {
        assert(getMaxBinNode(end - begin) <= binNode->END - binNode->BEGIN);
        if (end - begin <= n_leaf_limit) {
            if (direction) {
                for (S32 i = begin; i < end; i++) val_buf[i] = val[i];
            }
            *node = makeTreeNode(true, begin, end);
            return;
        }
        *node = makeTreeNode(false, begin, end);
        countBodies(val, begin, end, level, binNode);

        S32 offset[N_TREE_CHILDREN];
        exclusiveScan(offset, binNode->NBODY, begin);

        moveBodies(val, val_buf, begin, end, level, binNode, offset);

        BinaryTreeNode * binNodeOffset = binNode->BEGIN;
        BinaryTreeNode binNodeChild[N_TREE_CHILDREN];

        mk_task_group_w(end - begin);
        for (S32 i = 0; i < N_TREE_CHILDREN; i++) {
            S32 maxBinNode = getMaxBinNode(binNode->NBODY[i]);
            assert(binNodeOffset + maxBinNode <= binNode->END);
            binNodeChild[i].BEGIN = binNodeOffset;
            binNodeChild[i].END   = binNodeOffset + maxBinNode;
            S32 nbody = binNode->NBODY[i];
            auto lambda = [&, i, nbody] {
                sortBodiesImpl(val_buf, val, offset[i], offset[i] + nbody,
                               level + 1, &binNodeChild[i], n_leaf_limit,
                               !direction, &((*node)->CHILD[i]));
            };
            create_taskc_w(lambda, nbody);
            binNodeOffset += maxBinNode;
        }
        wait_tasks;
        for (S32 i = 0; i < N_TREE_CHILDREN; i++) {
            (*node)->NNODE += (*node)->CHILD[i]->NNODE;
        }
    }

    template<class Tobj>
    void sortBodies(Tobj * val, Tobj * val_buf, const S32 begin, const S32 end,
                    const S32 n_leaf_limit, TreeNode ** N0) {
        BinaryTreeNode binNode[1];
        int maxBinNode = getMaxBinNode(end - begin);
        binNode->BEGIN = new BinaryTreeNode[maxBinNode];
        binNode->END   = binNode->BEGIN + maxBinNode;
        sortBodiesImpl(val, val_buf, begin, end, 1, binNode, n_leaf_limit, false, N0);
        delete[] binNode->BEGIN;
    }

    template<class Ttc>
    void nodes2cells(TreeNode * node, Ttc * tc_first, const S32 adr_tc,
                     const S32 adr_child_tc, S32 * numLevels, const S32 level) {
        Ttc * cur_tc = tc_first + adr_tc;
        cur_tc->n_ptcl_   = node->NBODY;
        cur_tc->adr_ptcl_ = node->IBODY;
        cur_tc->level_    = level;
        cur_tc->adr_tc_   = adr_child_tc;
        if (node->NNODE <= 1) {
            *numLevels = std::max(*numLevels, level);
        } else {
            S32 adr_next_child_tc = adr_child_tc + N_TREE_CHILDREN;
            mk_task_group_w(node->NNODE);
            for (S32 i = 0; i < N_TREE_CHILDREN; i++) {
                auto lambda = [&, i, adr_next_child_tc] {
                    nodes2cells(node->CHILD[i], tc_first, adr_child_tc + i,
                                adr_next_child_tc, numLevels, level + 1);
                };
                create_taskc_w(lambda, node->CHILD[i]->NNODE);
                adr_next_child_tc += node->CHILD[i]->NNODE - 1;
            }
            wait_tasks;
            for (S32 i = 0; i < N_TREE_CHILDREN; i++) {
                delete node->CHILD[i];
            }
            *numLevels = std::max(*numLevels, level + 1);
        }
    }
#endif

    /*
    template<class T>
    class GetKeyForSort{
    public:
        U64 operator () (const T & val) const {
            return val.getKey();
        }
    };
    */

    // T is unsigned intger only
    template<class T, int NBIT=8>
    class RadixSort{
    private:
        int n_thread_;
        int n_bucket_;
        T mask_;
        int ** bucket_size_;
        int ** prefix_sum_;
    public:
        RadixSort(){
            static const unsigned long int one = 1;
#if defined(PARTICLE_SIMULATOR_THREAD_PARALLEL) && defined(_OPENMP)
            n_thread_ = omp_get_max_threads();
#else
            n_thread_ = 1;
#endif
            //std::cerr<<"n_thread_="<<n_thread_<<std::endl;
            n_bucket_ = 1<<NBIT;
            //std::cerr<<"n_bucket_="<<n_bucket_<<std::endl;
            mask_ = 0;
            for(int i=0; i<NBIT; i++) mask_ |= one<<i;
            //std::cerr<<std::hex<<"mask_="<<mask_<<std::endl;
            //std::cerr<<std::dec;
            bucket_size_ = new int * [n_thread_];
            prefix_sum_ = new int * [n_thread_];
            for(int ith=0; ith<n_thread_; ith++){
                bucket_size_[ith] = new int [n_bucket_];
                prefix_sum_[ith] = new int [n_bucket_];
            }
        }
        ~RadixSort(){
            for(int ith=0; ith<n_thread_; ith++){
                delete [] bucket_size_[ith];
                delete [] prefix_sum_[ith];
            }
            delete [] prefix_sum_;
            delete [] bucket_size_;
        }

        template<class Tobj>
        void lsdSort(Tobj * val,
                     Tobj * val_buf,
                     const int first,
                     const int last){
            //GetKeyForSort<Tobj> gk;
            //Tgk gk;
            const int n_tot = last - first + 1;
            int n_loop = (sizeof(T) * 8) / NBIT;
            if((sizeof(T) * 8) % NBIT != 0) n_loop++;

            for(int loop=0, shift=0; ; loop++, shift += NBIT){
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL		
#pragma omp parallel
#endif		
                {
#if defined(PARTICLE_SIMULATOR_THREAD_PARALLEL) && defined(_OPENMP)		    
                    int id_th = omp_get_thread_num();
#else
                    int id_th = 0;
#endif
                    for(int ibkt=0; ibkt<n_bucket_; ibkt++){
                        bucket_size_[id_th][ibkt] = 0;
                    }
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL		    
#pragma omp for
#endif
                    for(int i=0; i<n_tot; i++){
                        T radix = (val[i].getKey() >> shift) & mask_;
                        //T radix = ( gk(val[i]) >> shift) & mask_;
                        bucket_size_[id_th][radix]++;
                    }
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp single
#endif
                    {
                        prefix_sum_[0][0] = 0;
                        for(int j=1; j<n_thread_; j++){
                            prefix_sum_[j][0] = prefix_sum_[j-1][0] + bucket_size_[j-1][0];
                        }
                        for(int i=1; i<n_bucket_; i++){
                            prefix_sum_[0][i] = prefix_sum_[n_thread_-1][i-1] + bucket_size_[n_thread_-1][i-1];
                            for(int j=1; j<n_thread_; j++){
                                prefix_sum_[j][i] = prefix_sum_[j-1][i] + bucket_size_[j-1][i];
                            }
                        }
                    }
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL		    
#pragma omp for
#endif
                    for(int i=0; i<n_tot; i++){
                        T radix = (val[i].getKey() >> shift) & mask_;
                        //T radix = ( gk(val[i]) >> shift) & mask_;
                        const int offset = prefix_sum_[id_th][radix];
                        val_buf[offset] = val[i];
                        prefix_sum_[id_th][radix]++;
                    }
                } // OMP parallel scope
		if( loop+1 == n_loop) break;
                Tobj * Ptmp = val;
                val = val_buf;
                val_buf = Ptmp;
            }
            if(n_loop % 2 == 1){
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL		
#pragma omp parallel for
#endif
                for(int i=0; i<n_tot; i++){
                    val[i] = val_buf[i];
                }
            }
        }
    };
}

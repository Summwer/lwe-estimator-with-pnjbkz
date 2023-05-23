#include "enumbs.h"




struct bsNode
{
    EnumBS::blocksize_strategy bs;
    vector<bsNode> childs;
};


class SearchTree
{
public:
    // SearchTree(){root = Create(root);}
    EnumBS* enumbs;
    BKZJSim* sim;
    COST* cost;
    Params* params;
    SearchTree(Params* params,int d){
        sim = new BKZJSim(params,d);
        cost = new COST(params);
        enumbs = new EnumBS(params);
        this->params = params;
    }

    void search_tree_est(vector<double> l0);
    void search_tree_est_in_parallel(vector<double> l0);

    // ~SearchTree(){Release(root);}
    // void PreOrder(){PreOrder(root);}	//前序遍历
    // void InOrder(){InOrder(root);}		//中序遍历
    // void PostOrder(){PostOrder(root);}	//后序遍历

// private:
//     bsNode<blocksize_strategy> * root;
//     bsNode<blocksize_strategy> * Create(bsNode<blocksize_strategy> *bt);
//     void Release(bsNode<blocksize_strategy> *bt);
    
    // void PreOrder(bsNode<blocksize_strategy> *bt);
    // void InOrder(bsNode<blocksize_strategy> *bt);
    // void PostOrder(bsNode<blocksize_strategy> *bt);
};

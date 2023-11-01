#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <queue>

using namespace std;

int dimUsed;

struct KdNode
{
    KdNode *parent;
    KdNode *left;
    KdNode *right;
    std::vector<double> point;    
    int axis;       
    KdNode( std::vector<double> &data, const int &ax )
    {
        point = data;
        axis = ax;
        parent = nullptr;
        left = nullptr;
        right = nullptr;
    }
};

struct CompareKdNode {
    double distance_;
    KdNode* node;

    CompareKdNode( double dist, KdNode* n ) : distance_( dist ), node(n) {}

    bool operator>( const CompareKdNode& other ) const {
        return distance_ > other.distance_;
    }
};

double distance( const std::vector<double>& p1, const std::vector<double>& p2 ) {
    double dx = p1[0] - p2[0];
    double dy = p1[1] - p2[1];
    return sqrt( dx * dx + dy * dy );
}


ostream & operator<<(ostream & os, vector<double> vi)
{
    os << "(";
    for (int i = 0; i < vi.size(); i++)
        cout << vi[i] << ",";
    os << ")";
    return os;
}

class KdTree
{
private:
    std::vector<std::vector<double>> points_;
    KdNode* root_;
    int dimension_;

public:
    KdTree( const std::vector<std::vector<double>> &points, const int &dim ) 
        : points_( points ), dimension_( dim ) {}

    ~KdTree() {}

    void createTree()
    {
        root_ = createTreeNode( 0, points_.size() - 1, 0 );
    }

    KdNode* createTreeNode( const int &left, const int &right, const int &dim )
    {
        dimUsed = dim;
        if ( right < left ) return nullptr;
        // 按照k维进行排序
        
        sort( points_.begin() + left, points_.begin() + right + 1, [dim] ( const auto &p1, const auto &p2 ) { return p1[dim] < p2[dim]; });
        int mid = ( left + right + 1 ) / 2;
        KdNode* root = new KdNode( points_[mid], dim );
        root->left = createTreeNode( left, mid - 1, ( dim + 1 ) % dimension_);
        if ( root->left != nullptr )
            root->left->parent = root;
        root->right = createTreeNode( mid + 1, right, (dim + 1) % dimension_);
        if ( root->right != nullptr )
            root->right->parent = root;
        return root;
    }

    KdNode* searchKdTree( std::vector<double> &searchpoint ) 
    {
        int dim = 0;
        double minDis = numeric_limits<double>::max();
        KdNode* root = root_;
        KdNode* tmp;

        while ( root != nullptr )
        {
            tmp = root;
            if ( searchpoint[dim] < root->point[dim] )                           
                root = root->left;   
            else
                root = root->right;
            dim = ( dim + 1 ) % dimension_;
        }
        // 找到属于的那个节点
        cout<<"found"<<tmp->point << endl;
        minDis = min( distance( searchpoint, tmp->point ), minDis );
        KdNode* nearNode = tmp;
        // 退回到根节点
        while ( tmp->parent != nullptr )
        {
            tmp = tmp->parent;
            // 判断父节点是否更近，如果近，更新最近节点
            if ( distance( tmp->point, searchpoint ) < minDis )
            {
                nearNode = tmp;
                minDis = distance( tmp->point, searchpoint );
            }
            KdNode* son;
            // 判断当前轴与点的距离，如果小于minDis，则进行到另一半进行查找
            if ( abs( tmp->point[tmp->axis] - searchpoint[tmp->axis]) < minDis )
            {
                if( tmp->point[tmp->axis] > searchpoint[tmp->axis] )
                    son = tmp->right;
                else
                    son = tmp->left;
                searchKdTreeNode( searchpoint, minDis, nearNode, son );
            }
        }   
        // 对根节点的另外半边子树进行搜索
        /*if (abs(tmp->val[tmp->axis] - d[tmp->axis]) < minDis)
        {
            if (tmp->val[tmp->axis] > d[tmp->axis])
                tmp = tmp->rightChild;
            else
                tmp = tmp->leftChild;
            searchKdTreeNode(d, minDis, nearNode, tmp);
        }*/
        return nearNode;
    }

    void searchKdTreeNode( const std::vector<double> &searchpoint, double &minDis, KdNode* &nearNode, KdNode* tmp ) 
    {
    // 递归终止
        if ( tmp == nullptr ) return;
        // 判断当前节点是否小于
        if ( distance( tmp->point, searchpoint ) < minDis ) {
            minDis = distance( tmp->point, searchpoint );
            nearNode = tmp;
        }
        // 如果轴与节点的距离小于minDis，则两个半边都需要搜索，否则只需要搜索一个半边
        if ( abs( tmp->point[tmp->axis] - searchpoint[tmp->axis] ) < minDis ) {
            searchKdTreeNode( searchpoint, minDis, nearNode, tmp->left );
            searchKdTreeNode( searchpoint, minDis, nearNode, tmp->right );
        } else {
            // 选择搜索的一个半边
            if ( tmp->point[tmp->axis] > searchpoint[tmp->axis])
                searchKdTreeNode( searchpoint, minDis, nearNode, tmp->left );
            else
                searchKdTreeNode( searchpoint, minDis, nearNode, tmp->right );
        }
    }

    void printKdTree()
    {
        printKdTreeNode(root_);
    }

    void printKdTreeNode( KdNode * r )
    {
        if ( r == nullptr)
            return;
        printKdTreeNode(r->left);
        cout << r->point << "\t";
        printKdTreeNode(r->right);
    }

    // k-Nearest Neighbors search function
    std::vector<KdNode*> kNearestNeighbors(std::vector<double>& searchpoint, int k) {
        // Create a priority queue to store k-nearest neighbors
        std::priority_queue<pair<double, KdNode*>> nearestNeighbors;
        // std::priority_queue<CompareKdNode*> nearestNeighbors;

        // Perform k-NN search in the Kd-tree
        kNearestNeighborsHelper(root_, searchpoint, nearestNeighbors, k);

        // Extract the k-nearest neighbors from the priority queue
        std::vector<KdNode*> result;
        while (!nearestNeighbors.empty()) {
            result.push_back(nearestNeighbors.top().second);
            cout<<"einmal"<<endl;
            nearestNeighbors.pop();
        }
        return result;
    }

    void kNearestNeighborsHelper(KdNode* node, const std::vector<double>& searchpoint, std::priority_queue<pair<double, KdNode*>>& nearestNeighbors, int k) {
        if (node == nullptr)
            return;

        double distanceToNode = distance(searchpoint, node->point);

        if (nearestNeighbors.size() < k) {
            nearestNeighbors.push({distanceToNode, node});
        } else {
            if (distanceToNode < nearestNeighbors.top().first) {
                nearestNeighbors.pop();
                nearestNeighbors.push({distanceToNode, node});
            }
        }

        int axis = node->axis;
        double axisDifference = searchpoint[axis] - node->point[axis];

        if (axisDifference <= 0) {
            kNearestNeighborsHelper(node->left, searchpoint, nearestNeighbors, k);
            if (axisDifference * axisDifference < nearestNeighbors.top().first) {
                kNearestNeighborsHelper(node->right, searchpoint, nearestNeighbors, k);
            }
        } else {
            kNearestNeighborsHelper(node->right, searchpoint, nearestNeighbors, k);
            if (axisDifference * axisDifference < nearestNeighbors.top().first) {
                kNearestNeighborsHelper(node->left, searchpoint, nearestNeighbors, k);
            }
        }
    }
};


int main(int argc, char const *argv[])
{
    std::vector<std::vector<double>> points = { {2.0, 3.2}, {1.5, 12.1}, {1.15, 4.9}, {0.5, 3.9}, {2.3, 4.9}, {10.5, 4.3}, {7.5, 1.9}, {0.3, 2.7}, {0.9, 7.3} };
    KdTree* kdtree = new KdTree( points, 2);
    kdtree->createTree();
    kdtree->printKdTree();
    cout << endl;
    vector<double> vi = {2.0, 3.2};
    KdNode* r = kdtree->searchKdTree(vi);
    cout << r->point << endl;
    cout<<"----------------------------------------------------"<<endl;

    std::vector<double> searchPoint = {2.0, 3.2};

    // Perform k-NN search for k=3
    int k = 4;
    std::vector<KdNode*> nearestNeighbors = kdtree->kNearestNeighbors(searchPoint, k);

    // Print the k-NN results
    cout << "K-Nearest Neighbors for " << searchPoint << " with k=" << k << ":" << endl;
    for (int i = 0; i < k; i++) {
        cout << "Neighbor " << i + 1 << ": " << nearestNeighbors[i]->point << endl;
    }
    // cout<< nearestNeighbors.size();

    return 0;
}

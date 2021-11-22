#include <bits/stdc++.h>
using namespace std;
struct edge
{
    int first, second, weight;
};

class Graph
{
protected:
    int numberOfNodes, numberOfEdges;
    vector<vector<int>> listOfNeighbours;
    bool isDirected;

public:
    Graph(int numberOfNodes = 0, int numberOfEdges = 0, bool isDirected = 0, vector<vector<int>> listOfNeighbours = {})
    {
        this->numberOfNodes = numberOfNodes;
        this->numberOfEdges = numberOfEdges;
        this->isDirected = isDirected;
        this->listOfNeighbours = listOfNeighbours;
    }
    ~Graph()
    {
        this->numberOfEdges = 0;
        this->numberOfNodes = 0;
    }
    void setEdge(int, int);
    void printListOfNeighbours() const;
    int numberOfConnectedComponents(vector<int> &);
    vector<int> BFS(int, vector<int> &);
    stack<int> topologicalSorting(vector<int> &, stack<int> &);
    vector<vector<int>> stronglyConnectedComponents();
    vector<vector<int>> biconnectedComponents();
    bool HavelHakimi(vector<int>);
    vector<vector<int>> criticalConnections();
    vector<int> disjoint(vector<pair<int, pair<int, int>>>);
    int find(int, vector<int> &);
    void unite(int, int, vector<int> &);

private:
    void DFS(int, vector<int> &);
    void topologicalSortingDFS(int, vector<int> &, stack<int> &);
    void criticalDFS(int, int &, vector<int> &, vector<int> &, vector<int> &, int &, vector<vector<int>> &);
    void stronglyConnectedComponentsDFS(int, int &, vector<int> &, vector<int> &, vector<int> &, stack<int> &, vector<vector<int>> &);
    void biconnectedComponentsDFS(int, int &, int &, vector<int> &, vector<int> &, stack<int> &, vector<vector<int>> &);
};

class WeightedGraph : public Graph
{
private:
    vector<vector<int>> weightedEdges;
    vector<vector<edge>> listOfNeighboursWeight;

public:
    WeightedGraph(int numberOfNodes = 0, int numberOfEdges = 0, bool isDirected = 0, vector<vector<int>> listOfNeighbours = {}, vector<vector<int>> weightedEdges = {}) : Graph(numberOfNodes, numberOfEdges, isDirected, listOfNeighbours)
    {
        this->weightedEdges = weightedEdges;
    }
    void setWeightedEdge(int, int, int);
    vector<vector<int>> apm(vector<int> &, int &);
    vector<int> BellmanFord(int node);
    vector<int> Dijkstra(int node);
};
void WeightedGraph::setWeightedEdge(int firstNode, int secondNode, int weight)
{
    struct edge edgeValues;
    listOfNeighboursWeight.resize(this->numberOfNodes + 1);
    edgeValues.second = secondNode;
    edgeValues.weight = weight;
    listOfNeighboursWeight[firstNode].push_back(edgeValues);
    if (!this->isDirected)
        listOfNeighboursWeight[secondNode].push_back(edgeValues);
}

void Graph::setEdge(int firstNode, int secondNode)
{
    listOfNeighbours.resize(numberOfNodes + 1);
    listOfNeighbours[firstNode].push_back(secondNode);
    if (!isDirected)
        listOfNeighbours[secondNode].push_back(firstNode);
}
void Graph::printListOfNeighbours() const
{
    for (int i = 1; i < this->listOfNeighbours.size() - 1; i++)
    {
        cout << i << ": ";
        for (int j = 0; j < listOfNeighbours[i].size(); j++)
            cout << this->listOfNeighbours[i][j] << " ";
        cout << '\n';
    }
}
void Graph::DFS(int node, vector<int> &visited)
{
    visited[node] = 1;
    for (int neighbour : listOfNeighbours[node])
    {
        if (!visited[neighbour])
            DFS(neighbour, visited);
    }
}
int Graph::numberOfConnectedComponents(vector<int> &visited)
{
    int numberOfConnectedComponents = 0;
    for (int i = 1; i <= this->numberOfNodes; i++)
    {
        if (!visited[i])
        {
            numberOfConnectedComponents++;
            DFS(i, visited);
        }
    }
    return numberOfConnectedComponents;
}

vector<int> Graph::BFS(int sourceNode, vector<int> &visited)
{
    vector<int> distance(numberOfNodes + 1, 0);
    queue<int> queueBFS;
    queueBFS.push(sourceNode);
    visited[sourceNode] = 1;
    while (!queueBFS.empty())
    {
        for (int node : listOfNeighbours[queueBFS.front()])
        {
            if (!visited[node])
            {
                visited[node] = 1;
                queueBFS.push(node);
                distance[node] = distance[queueBFS.front()] + 1;
            }
        }
        queueBFS.pop();
    }
    return distance;
}

void Graph::topologicalSortingDFS(int node, vector<int> &visited, stack<int> &stack)
{
    visited[node] = 1;
    for (int neighbour : listOfNeighbours[node])
    {
        if (!visited[neighbour])
            topologicalSortingDFS(neighbour, visited, stack);
    }
    stack.push(node);
}

stack<int> Graph::topologicalSorting(vector<int> &visited, stack<int> &stack)
{
    for (int i = 1; i <= this->numberOfNodes; i++)
    {
        if (!visited[i])
            topologicalSortingDFS(i, visited, stack);
    }
    return stack;
}
void Graph::stronglyConnectedComponentsDFS(int node, int &currentLevel, vector<int> &level, vector<int> &lowestLevel, vector<int> &onStack, stack<int> &stack, vector<vector<int>> &solution)
{
    stack.push(node);
    onStack[node] = 1;
    level[node] = lowestLevel[node] = currentLevel;
    currentLevel++;

    for (auto neighbour : listOfNeighbours[node])
    {
        if (!level[neighbour])
            stronglyConnectedComponentsDFS(neighbour, currentLevel, level, lowestLevel, onStack, stack, solution);
        if (onStack[neighbour])
            lowestLevel[node] = min(lowestLevel[node], lowestLevel[neighbour]);
    }

    if (level[node] == lowestLevel[node])
    {
        vector<int> SCC;
        while (!stack.empty())
        {
            int top = stack.top();
            stack.pop();
            onStack[top] = 0;
            SCC.push_back(top);
            lowestLevel[top] = level[node];
            if (top == node)
                break;
        }
        solution.push_back(SCC);
    }
}

vector<vector<int>> Graph::stronglyConnectedComponents()
{

    int currentLevel = 1;
    stack<int> stack;
    vector<vector<int>> solution;
    vector<int> level(numberOfNodes + 1, 0);
    vector<int> lowestLevel(numberOfNodes + 1, 0);
    vector<int> onStack(numberOfNodes + 1, 0);

    for (int node = 1; node <= numberOfNodes; node++)
        if (level[node] == 0)
            stronglyConnectedComponentsDFS(node, currentLevel, level, lowestLevel, onStack, stack, solution);

    return solution;
}
void Graph::biconnectedComponentsDFS(int node, int &parent, int &currentLevel, vector<int> &level, vector<int> &lowestLevel, stack<int> &stack, vector<vector<int>> &solution)
{
    stack.push(node);
    level[node] = lowestLevel[node] = currentLevel;
    currentLevel++;

    for (auto neighbour : listOfNeighbours[node])
    {
        if (neighbour == parent)
            continue;

        else if (!level[neighbour])
        {

            biconnectedComponentsDFS(neighbour, node, currentLevel, level, lowestLevel, stack, solution);
            lowestLevel[node] = min(lowestLevel[node], lowestLevel[neighbour]);
            if (lowestLevel[neighbour] >= level[node])
            {
                vector<int> BCC;
                while (stack.top() != neighbour)
                {
                    BCC.push_back(stack.top());
                    stack.pop();
                }
                BCC.push_back(neighbour);
                stack.pop();
                BCC.push_back(node);
                solution.push_back(BCC);
            }
        }
        else
            lowestLevel[node] = min(lowestLevel[node], level[neighbour]);
    }
}
vector<vector<int>> Graph::biconnectedComponents()
{
    int currentLevel = 1;
    stack<int> stack;
    vector<vector<int>> solution;
    vector<int> level(numberOfNodes + 1, 0);
    vector<int> lowestLevel(numberOfNodes + 1, 0);
    int parent = 0;
    for (int node = 1; node <= numberOfNodes; node++)
        if (level[node] == 0)
            biconnectedComponentsDFS(node, parent, currentLevel, level, lowestLevel, stack, solution);

    return solution;
}
bool Graph::HavelHakimi(vector<int> degrees)
{

    sort(degrees.begin(), degrees.end(), greater<>());

    while (1)
    {
        int node = degrees[0];
        if (node == 0)
            return 1;
        degrees.erase(degrees.begin() + 0);

        if (degrees.size() < node)
            return 0;

        for (int i = 0; i < node; i++)
        {
            degrees[i]--;

            if (degrees[i] < 0)
                return 0;
        }
    }
}
void Graph::criticalDFS(int node, int &parent, vector<int> &visited, vector<int> &lowestLevel, vector<int> &level, int &currentLevel, vector<vector<int>> &criticalEdges)
{
    visited[node] = 1;
    lowestLevel[node] = level[node] = currentLevel;
    currentLevel++;
    for (auto neighbour : listOfNeighbours[node])
    {
        if (neighbour == parent)
            continue;

        if (visited[neighbour] == 1)
            lowestLevel[node] = min(lowestLevel[node], level[neighbour]);
        else
        {
            criticalDFS(neighbour, node, visited, lowestLevel, level, currentLevel, criticalEdges);
            if (lowestLevel[neighbour] > level[node])
                criticalEdges.push_back({node, neighbour});

            lowestLevel[node] = min(lowestLevel[node], lowestLevel[neighbour]);
        }
    }
}
vector<vector<int>> Graph::criticalConnections()
{
    vector<int> level(numberOfNodes + 1, 0);
    vector<int> lowestLevel(numberOfNodes + 1, 0);
    vector<int> visited(numberOfNodes + 1, 0);
    vector<vector<int>> criticalEdges;
    int currentLevel, parent;
    currentLevel = 1;
    parent = 0;
    criticalDFS(1, parent, visited, lowestLevel, level, currentLevel, criticalEdges);
    return criticalEdges;
}
int Graph::find(int node, vector<int> &root)
{
    if (root[node] != node)
        root[node] = find(root[node], root);

    return root[node];
}
void Graph::unite(int x, int y, vector<int> &root)
{
    int firstRoot, secondRoot;
    firstRoot = find(x, root);
    secondRoot = find(y, root);
    root[firstRoot] = secondRoot;
}

vector<int> Graph::disjoint(vector<pair<int, pair<int, int>>> tasks)
{
    vector<int> root(numberOfNodes + 1);
    vector<int> solution;
    for (int i = 1; i <= numberOfNodes; i++)
        root[i] = i;
    for (int i = 0; i < tasks.size(); i++)
    {
        int op = tasks[i].first;
        int x = tasks[i].second.first;
        int y = tasks[i].second.second;

        if (op == 1)
            unite(x, y, root);

        else if (op == 2)
        {
            if (find(x, root) == find(y, root))
                solution.push_back(1);
            else
                solution.push_back(0);
        }
    }
    return solution;
}

bool comp(vector<int> first, vector<int> second)
{
    return first[2] < second[2];
}
vector<vector<int>> WeightedGraph::apm(vector<int> &parent, int &totalWeight)
{
    vector<vector<int>> solution;
    int firstNode, secondNode, firstParent, secondParent, i = 0;
    sort(weightedEdges.begin(), weightedEdges.end(), comp);
    while (solution.size() < numberOfNodes - 1)
    {
        firstNode = weightedEdges[i][0];
        secondNode = weightedEdges[i][1];
        firstParent = firstNode;
        secondParent = secondNode;

        firstParent = find(firstNode, parent);
        secondParent = find(secondNode, parent);
        if (firstParent != secondParent)
        {
            if (firstParent < secondParent)
                parent[firstParent] = parent[secondParent];

            else
                parent[secondParent] = parent[firstParent];

            solution.push_back({firstNode, secondNode});
            totalWeight += weightedEdges[i][2];
        }
        i++;
    }
    return solution;
}
vector<int> WeightedGraph::BellmanFord(int node)
{
    vector<int> inQueue(numberOfNodes + 1, 0);
    vector<int> minWeight(numberOfNodes + 1, INT_MAX);
    vector<int> visited(numberOfNodes + 1, 0);
    queue<int> queue;

    minWeight[node] = 0;
    queue.push(node);
    inQueue[node] = 1;

    while (!queue.empty())
    {
        int currentNode = queue.front();
        visited[currentNode]++;

        if (visited[currentNode] >= numberOfNodes)
            return {-1};
        queue.pop();
        inQueue[currentNode] = 0;

        for (auto it : listOfNeighboursWeight[currentNode])
        {
            int neighbour = it.second;
            int weight = it.weight;

            if (minWeight[currentNode] + weight < minWeight[neighbour])
            {
                minWeight[neighbour] = minWeight[currentNode] + weight;
                if (!inQueue[neighbour])
                {
                    inQueue[neighbour] = 1;
                    queue.push(neighbour);
                }
            }
        }
    }
    return minWeight;
}
vector<int> WeightedGraph::Dijkstra(int node)
{
    vector<int> minWeight(numberOfNodes + 1, INT_MAX);
    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> pQueue;

    minWeight[node] = 0;
    pQueue.push({node, 0});

    while (!pQueue.empty())
    {
        int currentNode = pQueue.top().first;
        pQueue.pop();

        for (auto it : listOfNeighboursWeight[currentNode])
        {
            int neighbour = it.second;
            int weight = it.weight;

            if (minWeight[currentNode] + weight < minWeight[neighbour])
            {
                minWeight[neighbour] = minWeight[currentNode] + weight;
                pQueue.push({neighbour, minWeight[neighbour]});
            }
        }
    }
    for (int i = 2; i <= numberOfNodes; i++)
        if (minWeight[i] == INT_MAX)
            minWeight[i] = 0;

    return minWeight;
}
void DFSinfoarena()
{
    ifstream in("dfs.in");
    ofstream out("dfs.out");
    int numberOfNodes, numberOfEdges, firstNode, secondNode;
    in >> numberOfNodes >> numberOfEdges;
    Graph Gr(numberOfNodes, numberOfEdges, 0);
    for (int i = 1; i <= numberOfEdges; i++)
    {
        in >> firstNode >> secondNode;
        Gr.setEdge(firstNode, secondNode);
    }
    vector<int> visited(numberOfNodes + 1, 0);
    out << Gr.numberOfConnectedComponents(visited);
}
void BFSinfoarena()
{
    ifstream in("bfs.in");
    ofstream out("bfs.out");
    int numberOfNodes, numberOfEdges, firstNode, secondNode, sourceNode;
    in >> numberOfNodes >> numberOfEdges >> sourceNode;
    Graph Gr(numberOfNodes, numberOfEdges, 1);
    for (int i = 1; i <= numberOfEdges; i++)
    {
        in >> firstNode >> secondNode;
        Gr.setEdge(firstNode, secondNode);
    }
    vector<int> visited(numberOfNodes + 1, 0);
    vector<int> distance;
    distance = Gr.BFS(sourceNode, visited);

    for (int i = 1; i < distance.size(); i++)
    {
        if (distance[i] == 0 && sourceNode != i)
            out << "-1 ";
        else
            out << distance[i] << " ";
    }
}
void sortaretInfoarena()
{
    ifstream in("sortaret.in");
    ofstream out("sortaret.out");
    int numberOfNodes, numberOfEdges, firstNode, secondNode;
    in >> numberOfNodes >> numberOfEdges;
    Graph Gr(numberOfNodes, numberOfEdges, 1);
    for (int i = 1; i <= numberOfEdges; i++)
    {
        in >> firstNode >> secondNode;
        Gr.setEdge(firstNode, secondNode);
    }
    vector<int> visited(numberOfNodes + 1, 0);
    stack<int> stack;
    stack = Gr.topologicalSorting(visited, stack);
    while (!stack.empty())
    {
        out << stack.top() << " ";
        stack.pop();
    }
}
void ctcInfoarena()
{
    ifstream in("ctc.in");
    ofstream out("ctc.out");
    int numberOfNodes, numberOfEdges, firstNode, secondNode;
    in >> numberOfNodes >> numberOfEdges;
    Graph Gr(numberOfNodes, numberOfEdges, 1);
    for (int i = 1; i <= numberOfEdges; i++)
    {
        in >> firstNode >> secondNode;
        Gr.setEdge(firstNode, secondNode);
    }
    vector<vector<int>> solution = Gr.stronglyConnectedComponents();
    out << solution.size() << '\n';
    for (int i = 0; i < solution.size(); i++)
    {
        for (int j = 0; j < solution[i].size(); j++)
            out << solution[i][j] << " ";
        out << '\n';
    }
}
void biconexInfoarena()
{
    ifstream in("biconex.in");
    ofstream out("biconex.out");
    int numberOfNodes, numberOfEdges, firstNode, secondNode;
    in >> numberOfNodes >> numberOfEdges;
    Graph Gr(numberOfNodes, numberOfEdges, 0);
    for (int i = 1; i <= numberOfEdges; i++)
    {
        in >> firstNode >> secondNode;
        Gr.setEdge(firstNode, secondNode);
    }
    vector<vector<int>> solution = Gr.biconnectedComponents();
    out << solution.size() << '\n';
    for (int i = 0; i < solution.size(); i++)
    {
        for (int j = 0; j < solution[i].size(); j++)
            out << solution[i][j] << " ";
        out << '\n';
    }
}
void CCNleetcode()
{
    ifstream in("graph.in");
    ofstream out("graph.out");
    int numberOfNodes, numberOfEdges, firstNode, secondNode;
    in >> numberOfNodes;
    Graph Gr(numberOfNodes, numberOfNodes, 0);
    for (int i = 1; i <= numberOfNodes; i++)
    {
        in >> firstNode >> secondNode;
        Gr.setEdge(firstNode, secondNode);
    }
    vector<vector<int>> criticalEdges = Gr.criticalConnections();

    for (int i = 0; i < criticalEdges.size(); i++)
    {
        for (int j = 0; j < criticalEdges[i].size(); j++)
            out << criticalEdges[i][j] << " ";
        out << '\n';
    }
}
void HavelHakimiRez()
{
    ifstream in("graph.in");
    ofstream out("graph.out");
    int number, degree;
    vector<int> degrees(number + 1);
    for (int i = 0; i < number; i++)
    {
        in >> degree;
        degrees.push_back(degree);
    }

    //  out << Gr.HavelHakimi(degrees);
}
void disjointInfoarena()
{
    ifstream in("disjoint.in");
    ofstream out("disjoint.out");
    int numberOfSets, numberOfTasks, op, x, y;
    in >> numberOfSets >> numberOfTasks;
    Graph Gr(numberOfSets);
    vector<pair<int, pair<int, int>>> tasks;

    for (int i = 0; i < numberOfTasks; i++)
    {
        in >> op >> x >> y;
        tasks.push_back(make_pair(op, make_pair(x, y)));
    }
    vector<int> solution = Gr.disjoint(tasks);
    for (int i = 0; i < solution.size(); i++)
        if (solution[i])
            out << "DA" << '\n';
        else
            out << "NU" << '\n';
}
void apmInfoarena()
{
    ifstream in("apm.in");
    ofstream out("apm.out");
    int numberOfNodes, numberOfEdges, firstNode, secondNode, weight;
    in >> numberOfNodes >> numberOfEdges;
    vector<vector<int>> listOfNeighbours(numberOfNodes + 1, {0});
    vector<vector<int>> weightedEdges;
    for (int i = 0; i < numberOfEdges; i++)
    {
        in >> firstNode >> secondNode >> weight;
        weightedEdges.push_back({firstNode, secondNode, weight});
    }

    WeightedGraph Gr(numberOfNodes, numberOfEdges, 0, listOfNeighbours, weightedEdges);

    vector<int> parent(numberOfNodes + 1);
    for (int i = 1; i < parent.size(); i++)
        parent[i] = i;
    int totalWeight = 0;
    vector<vector<int>> solution = Gr.apm(parent, totalWeight);
    out << totalWeight << '\n'
        << solution.size() << '\n';
    for (int i = 0; i < solution.size(); i++)
        out << solution[i][0] << " " << solution[i][1] << '\n';
}
void BellmanFordInfoarena()
{
    ifstream in("bellmanford.in");
    ofstream out("bellmanford.out");
    int numberOfNodes, numberOfEdges, firstNode, secondNode, weight;
    in >> numberOfNodes >> numberOfEdges;

    WeightedGraph Gr(numberOfNodes, numberOfEdges, 1);
    for (int i = 0; i < numberOfEdges; i++)
    {
        in >> firstNode >> secondNode >> weight;
        Gr.setWeightedEdge(firstNode, secondNode, weight);
    }

    vector<int> minWeight = Gr.BellmanFord(1);
    if (minWeight[0] == -1)
        out << "Ciclu negativ!";
    else
        for (int i = 2; i < minWeight.size(); i++)
            out << minWeight[i] << " ";
}
void DijkstraInfoarena()
{
    ifstream in("dijkstra.in");
    ofstream out("dijkstra.out");
    int numberOfNodes, numberOfEdges, firstNode, secondNode, weight;
    in >> numberOfNodes >> numberOfEdges;

    WeightedGraph Gr(numberOfNodes, numberOfEdges, 1);
    for (int i = 0; i < numberOfEdges; i++)
    {
        in >> firstNode >> secondNode >> weight;
        Gr.setWeightedEdge(firstNode, secondNode, weight);
    }

    vector<int> minWeight = Gr.Dijkstra(1);
    for (int i = 2; i < minWeight.size(); i++)
        out << minWeight[i] << " ";
}
int main()
{
    //DFSinfoarena();
    //BFSinfoarena();
    //sortaretInfoarena();
    //ctcInfoarena();
    //biconexInfoarena();
    //CCNleetcode();
    //HavelHakimiRez();
    //disjointInfoarena();
    //BellmanFordInfoarena();
    //DijkstraInfoarena();
    return 0;
}

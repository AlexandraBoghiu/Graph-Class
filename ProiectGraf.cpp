#include <bits/stdc++.h>
using namespace std;
ifstream in("graph.in");
ofstream out("graph.out");

class Graph
{
private:
    int numberOfNodes, numberOfEdges;
    vector<vector<int>> listOfNeighbours;
    bool isDirected;

public:
    Graph(int numberOfNodes = 0, int numberOfEdges = 0, bool isDirected = 0)
    {
        this->numberOfNodes = numberOfNodes;
        this->numberOfEdges = numberOfEdges;
        this->isDirected = isDirected;
    }
    ~Graph()
    {
        this->numberOfEdges = 0;
        this->numberOfNodes = 0;
    }
    void setNumberOfEdges(int &);
    void setNumberOfNodes(int &);
    void setEdges(int, int);
    int getNumberOfNodes();
    int getNumberOfEdges();
    vector<int> getDegrees();
    void readEdges(istream &in);
    void printListOfNeighbours(ostream &out);
    int numberOfConnectedComponents(bool[]);
    void BFS(int, bool[], int[]);
    void topologicalSorting(bool[], stack<int> &);
    void stronglyConnectedComponents();
    void biconnectedComponents();
    bool HavelHakimi(vector<int>);
    void criticalConnections();

private:
    void DFS(int, bool[]);
    void topologicalSortingDFS(int, bool[], stack<int> &);
    void criticalDFS(int, int &, vector<int> &, vector<int> &, vector<int> &, int &, vector<vector<int>> &);
    void stronglyConnectedComponentsDFS(int, int &, vector<int> &, vector<int> &, bool[], stack<int> &, vector<vector<int>> &);
    void biconnectedComponentsDFS(int, int &, int &, vector<int> &, vector<int> &, stack<int> &, vector<vector<int>> &);
};

void Graph::setNumberOfEdges(int &m)
{
    numberOfEdges = m;
}
void Graph::setNumberOfNodes(int &n)
{
    numberOfNodes = n;
}
void Graph::setEdges(int firstNode, int secondNode)
{
    listOfNeighbours[firstNode].push_back(secondNode);
}
int Graph::getNumberOfNodes()
{
    return numberOfNodes;
}
int Graph::getNumberOfEdges()
{
    return numberOfEdges;
}
vector<int> Graph::getDegrees()
{
    vector<int> degrees(numberOfNodes + 1);
    for(int i = 0; i < this->listOfNeighbours.size(); i++)
        degrees[i] = listOfNeighbours[i].size();
    
    return degrees;
}
void Graph::readEdges(istream &in)
{
    int firstNode, secondNode;
    listOfNeighbours.resize(numberOfNodes + 2);
    for (int i = 0; i < this->numberOfEdges; i++)
    {
        in >> firstNode >> secondNode;
        setEdges(firstNode, secondNode);
        if (!isDirected)
            setEdges(secondNode, firstNode);
    }
}
void Graph::printListOfNeighbours(ostream &out)
{
    for (int i = 1; i < this->listOfNeighbours.size() - 1; i++)
    {
        out << i << ": ";
        for (int j = 0; j < listOfNeighbours[i].size(); j++)
            out << this->listOfNeighbours[i][j] << " ";
    }
}
void Graph::DFS(int node, bool visited[])
{
    visited[node] = 1;
    for (int neighbour : listOfNeighbours[node])
    {
        if (!visited[neighbour])
            DFS(neighbour, visited);
    }
}
int Graph::numberOfConnectedComponents(bool visited[])
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

void Graph::BFS(int sourceNode, bool visited[], int distance[])
{

    queue<int> queueBFS;
    queueBFS.push(sourceNode);
    visited[sourceNode] = 1;
    distance[numberOfNodes] = {0};
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
    for (int node = 1; node <= numberOfNodes; node++)
    {
        if (sourceNode != node && distance[node] == 0)
            out << "-1 ";
        else
            out << distance[node] << " ";
    }
}

void Graph::topologicalSortingDFS(int node, bool visited[], stack<int> &stack)
{
    visited[node] = 1;
    for (int neighbour : listOfNeighbours[node])
    {
        if (!visited[neighbour])
            topologicalSortingDFS(neighbour, visited, stack);
    }
    stack.push(node);
}

void Graph::topologicalSorting(bool visited[], stack<int> &stack)
{
    for (int i = 1; i <= this->numberOfNodes; i++)
    {
        if (!visited[i])
            topologicalSortingDFS(i, visited, stack);
    }

    while (!stack.empty())
    {
        out << stack.top() << " ";
        stack.pop();
    }
}
void Graph::stronglyConnectedComponentsDFS(int node, int &currentLevel, vector<int> &level, vector<int> &lowestLevel, bool onStack[], stack<int> &stack, vector<vector<int>> &solution)
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

void Graph::stronglyConnectedComponents()
{

    int currentLevel = 1;
    stack<int> stack;
    vector<vector<int>> solution;
    vector<int> level(numberOfNodes + 1, 0);
    vector<int> lowestLevel(numberOfNodes + 1, 0);
    bool onStack[numberOfNodes + 1] = {0};

    for (int node = 1; node <= numberOfNodes; node++)
        if (level[node] == 0)
            stronglyConnectedComponentsDFS(node, currentLevel, level, lowestLevel, onStack, stack, solution);

    out << solution.size() << '\n';
    for (int i = 0; i < solution.size(); i++)
    {
        for (int j = 0; j < solution[i].size(); j++)
            out << solution[i][j] << " ";
        out << '\n';
    }
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
void Graph::biconnectedComponents()
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

    out << solution.size() << '\n';
    for (int i = 0; i < solution.size(); i++)
    {
        for (int j = 0; j < solution[i].size(); j++)
            out << solution[i][j] << " ";
        out << '\n';
    }
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
void Graph::criticalConnections()
{
    vector<int> level(numberOfNodes + 1, 0);
    vector<int> lowestLevel(numberOfNodes + 1, 0);
    vector<int> visited(numberOfNodes + 1, 0);
    vector<vector<int>> criticalEdges;
    int currentLevel, parent;
    currentLevel = 1;
    parent = 0;
    criticalDFS(1, parent, visited, lowestLevel, level, currentLevel, criticalEdges);
    for (int i = 0; i < criticalEdges.size(); i++)
    {
        for (int j = 0; j < criticalEdges[i].size(); j++)
            out << criticalEdges[i][j] << " ";
        out << '\n';
    }
}
void DFSinfoarena()
{
    int numberOfNodes, numberOfEdges;
    in >> numberOfNodes >> numberOfEdges;
    Graph Gr(numberOfNodes, numberOfEdges, 0);
    Gr.readEdges(in);
    bool visited[numberOfNodes + 1] = {0};
    out << Gr.numberOfConnectedComponents(visited);
}
void BFSinfoarena()
{
    int numberOfNodes, numberOfEdges, sourceNode;
    in >> numberOfNodes >> numberOfEdges >> sourceNode;
    Graph Gr(numberOfNodes, numberOfEdges, 1);
    Gr.readEdges(in);
    bool visited[numberOfNodes + 1] = {0};
    int distance[numberOfNodes + 1] = {0};
    Gr.BFS(sourceNode, visited, distance);
}
void sortaretInfoarena()
{
    int numberOfNodes, numberOfEdges, sourceNode;
    in >> numberOfNodes >> numberOfEdges;
    Graph Gr(numberOfNodes, numberOfEdges, 1);
    Gr.readEdges(in);
    bool visited[numberOfNodes + 1] = {0};
    stack<int> stack;
    Gr.topologicalSorting(visited, stack);
}
void ctcInfoarena()
{
    int numberOfNodes, numberOfEdges;
    in >> numberOfNodes >> numberOfEdges;
    Graph Gr(numberOfNodes, numberOfEdges, 1);
    Gr.readEdges(in);
    Gr.stronglyConnectedComponents();
}
void biconexInfoarena()
{
    int numberOfNodes, numberOfEdges;
    in >> numberOfNodes >> numberOfEdges;
    Graph Gr(numberOfNodes, numberOfEdges, 0);
    Gr.readEdges(in);
    Gr.biconnectedComponents();
}
void CCNleetcode()
{
    int numberOfNodes, numberOfEdges;
    in >> numberOfNodes >> numberOfEdges;
    Graph Gr(numberOfNodes, numberOfEdges, 0);
    Gr.readEdges(in);
    Gr.criticalConnections();
}
void HavelHakimiRez()
{
    int numberOfNodes, numberOfEdges;
    vector<int> degrees(numberOfNodes + 1);
    in >> numberOfNodes >> numberOfEdges;
    Graph Gr(numberOfNodes, numberOfEdges, 0);
    Gr.readEdges(in);
    degrees = Gr.getDegrees();
    out<<Gr.HavelHakimi(degrees);
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
    return 0;
}

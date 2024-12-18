#include "7hw.h"


struct Edge {
    Node *node;
    double weight;

    Edge(Node *node, double weight) : node(node), weight(weight) {
    }
};

struct Node {
    double latitute;
    double longitute;
    std::vector<Edge> children;
    std::vector<Edge> parents;

    Node(const double latitute, const double longitute) : latitute(latitute), longitute(longitute) {
    }
};

struct Graph {
    std::vector<Node *> nodes;

    Node *get_node_ptr(const double latitute, const double longitute,
                       std::unordered_map<std::string, Node *> &nodes_map) {
        std::string key = std::to_string(latitute) + "," + std::to_string(longitute);
        if (nodes_map.find(key) != nodes_map.end())
            return nodes_map[key];

        Node *node_ptr = new Node(latitute, longitute);
        nodes.push_back(node_ptr);
        nodes_map[key] = node_ptr;
        return node_ptr;
    }

    std::pair<double, double> read_coords(const std::string str) {
        size_t comma = str.find(',');
        double longitute = std::stod(str.substr(0, comma));
        double latitute = std::stod(str.substr(comma + 1));
        return {latitute, longitute};
    }

    std::pair<std::pair<double, double>, double> read_edge(const std::string str) {
        size_t comma1 = str.find(',');
        size_t comma2 = str.find(',', comma1 + 1);
        double longitute = std::stod(str.substr(0, comma1));
        double latitute = std::stod(str.substr(comma1 + 1, comma2 - comma1 - 1));
        double weight = std::stod(str.substr(comma2 + 1));
        return {{latitute, longitute}, weight};
    }

    Graph(const std::string str) {
        std::ifstream file(str);
        std::string line;

        std::unordered_map<std::string, Node *> nodes_map;

        while (file >> line) {
            size_t first_sep = line.find(':');

            std::string node0_str = line.substr(0, first_sep);

            auto coords = read_coords(node0_str);
            Node *node0_ptr = get_node_ptr(coords.first, coords.second, nodes_map);

            while (line.find(';', first_sep + 1) != std::string::npos) {
                size_t second_sep = line.find(';', first_sep + 1);
                std::string node_str = line.substr(first_sep + 1, second_sep - first_sep - 1);
                auto edge1 = read_edge(node_str);
                Node *node1_ptr = get_node_ptr(edge1.first.first, edge1.first.second, nodes_map);
                node0_ptr->children.push_back({node1_ptr, edge1.second});
                node1_ptr->parents.push_back({node0_ptr, edge1.second});
                first_sep = second_sep;
            }
        }
    }
};

struct Pathway {
    std::vector<Node *> path;
    double distance;

    Pathway() : path(std::vector<Node *>()), distance(1e9) {
    }

    Pathway(std::vector<Node *> path, double distance) : path(path), distance(distance) {
    }
};

struct Algorithm_result {
    Pathway pathway;
    double time;

    Algorithm_result(std::vector<Node *> path, double distance, double time)
        : pathway(path, distance), time(time) {
    }
};

struct My_way {
    Node *start;
    Node *destination;
    Graph *graph;
    std::vector<Node *> extra_way_points;

    Node *closest_node(const double latitute, const double longitute) {
        Node *closest = nullptr;
        double min_distance = 1e9;
        for (auto node: graph->nodes) {
            double distance = (node->latitute - latitute) * (node->latitute - latitute) +
                              (node->longitute - longitute) * (node->longitute - longitute);
            if (distance < min_distance) {
                min_distance = distance;
                closest = node;
            }
        }
        return closest;
    }

    My_way(const double latitute_start, const double longitute_start,
           const double latitute_destination, const double longitute_destination,
           const std::string filename) {
        graph = new Graph(filename);
        start = closest_node(latitute_start, longitute_start);
        destination = closest_node(latitute_destination, longitute_destination);
    }

    void add_extra_point(const double latitute, const double longitute) {
        extra_way_points.push_back(closest_node(latitute, longitute));
    }

    std::vector<Node *> get_path(Node *s, Node *e, std::unordered_map<Node *, double> &dist) {
        std::vector<Node *> path;
        Node *cur = e;
        while (cur != s) {
            path.push_back(cur);
            for (auto edge: cur->parents)
                if (dist[cur] == dist[edge.node] + edge.weight) {
                    cur = edge.node;
                    break;
                }
        }
        path.push_back(s);
        std::reverse(path.begin(), path.end());
        return path;
    }

    Pathway dfs(Node *s, Node *e) {
        std::unordered_map<Node *, double> dist;
        for (auto node: graph->nodes)
            dist[node] = 1e9;
        std::stack<Node *> node_stack;
        node_stack.push(s);
        dist[s] = 0;

        while (!node_stack.empty()) {
            Node *node = node_stack.top();
            node_stack.pop();
            for (auto edge: node->children)
                if (dist[edge.node] > dist[node] + edge.weight) {
                    dist[edge.node] = dist[node] + edge.weight;
                    node_stack.push(edge.node);
                }
        }

        if (dist[e] == 1e9)
            return {};

        return {get_path(s, e, dist), dist[e]};
    }

    Pathway bfs(Node *s, Node *e) {
        std::unordered_map<Node *, double> dist;
        for (auto node: graph->nodes)
            dist[node] = 1e9;
        std::queue<Node *> node_queue;
        dist[s] = 0;
        node_queue.push(s);

        while (!node_queue.empty()) {
            Node *node = node_queue.front();
            node_queue.pop();
            for (auto edge: node->children)
                if (dist[edge.node] > dist[node] + edge.weight) {
                    dist[edge.node] = dist[node] + edge.weight;
                    node_queue.push(edge.node);
                }
        }

        if (dist[e] == 1e9)
            return {};

        return {get_path(s, e, dist), dist[e]};
    }

    Pathway dijkstra(Node *s, Node *e) {
        std::unordered_map<Node *, double> dist;
        for (auto node: graph->nodes)
            dist[node] = 1e9;
        auto cmp = [](std::pair<double, Node *> left, std::pair<double, Node *> right) {
            return left.first > right.first;
        };
        std::priority_queue<std::pair<double, Node *>, std::vector<std::pair<double, Node *> >, decltype(cmp)>
                node_queue(cmp);
        dist[s] = 0;
        node_queue.push({0, s});

        while (!node_queue.empty()) {
            Node *node = node_queue.top().second;
            node_queue.pop();
            if (node == e)
                break;
            for (auto edge: node->children)
                if (dist[edge.node] > dist[node] + edge.weight) {
                    dist[edge.node] = dist[node] + edge.weight;
                    node_queue.push({dist[edge.node], edge.node});
                }
        }

        if (dist[e] == 1e9)
            return {};

        return {get_path(s, e, dist), dist[e]};
    }

    double heuristic(Node *a, Node *b) {
        return (a->latitute - b->latitute) * (a->latitute - b->latitute) +
               (a->longitute - b->longitute) * (a->longitute - b->longitute);
    }

    Pathway a_star(Node *s, Node *e) {
        std::unordered_map<Node *, double> dist;
        for (auto node: graph->nodes)
            dist[node] = 1e9;
        auto cmp = [](std::pair<double, Node *> left, std::pair<double, Node *> right) {
            return left.first > right.first;
        };
        std::priority_queue<std::pair<double, Node *>, std::vector<std::pair<double, Node *> >, decltype(cmp)>
                node_queue(cmp);
        dist[s] = 0;
        node_queue.push({0, s});
        while (!node_queue.empty()) {
            Node *node = node_queue.top().second;
            node_queue.pop();
            if (node == e)
                break;
            for (auto edge: node->children)
                if (dist[edge.node] > dist[node] + edge.weight) {
                    dist[edge.node] = dist[node] + edge.weight;
                    node_queue.push({dist[edge.node] + heuristic(edge.node, e), edge.node});
                }
        }

        if (dist[e] == 1e9)
            return {};

        return {get_path(s, e, dist), dist[e]};
    }

    Pathway combinatory_pathway(std::set<size_t> points_left, size_t current_node,
                                std::vector<Pathway> &path_from_start,
                                std::vector<Pathway> &path_between_points, std::vector<Pathway> &path_to_destination) {
        if (points_left.empty())
            return path_to_destination[current_node - 1];

        Pathway best_pathway;
        size_t best_point = 0;

        if (current_node == 0) {
            for (auto point: points_left) {
                std::set<size_t> new_points_left = points_left;
                new_points_left.erase(point);

                Pathway current = combinatory_pathway(new_points_left, point,
                                                      path_from_start, path_between_points,
                                                      path_to_destination);
                if (current.distance + path_from_start[point - 1].distance <
                    best_pathway.distance + path_from_start[best_point - 1].distance) {
                    best_pathway = current;
                    best_point = point;
                }
            }
            Pathway result = path_from_start[best_point - 1];
            result.path.insert(result.path.end(), best_pathway.path.begin(), best_pathway.path.end());
            result.distance += best_pathway.distance;
            return result;
        }

        for (auto point: points_left) {
            std::set<size_t> new_points_left = points_left;
            new_points_left.erase(point);

            Pathway current = combinatory_pathway(new_points_left, point,
                                                  path_from_start, path_between_points,
                                                  path_to_destination);

            if (current.distance
                + path_between_points[(current_node - 1) * (extra_way_points.size() - 1) + point - 1].distance <
                best_pathway.distance
                + path_between_points[(current_node - 1) * (extra_way_points.size() - 1) + best_point - 1].distance) {
                best_pathway = current;
                best_point = point;
            }
        }
        Pathway result = path_between_points[(current_node - 1) * (extra_way_points.size() - 1) + best_point - 1];
        result.path.insert(result.path.end(), best_pathway.path.begin(), best_pathway.path.end());
        result.distance += best_pathway.distance;
        return result;
    }

    template<typename Function>
    Pathway find_best_path(Function algorithm) {
        if (extra_way_points.empty())
            return (this->*algorithm)(start, destination);

        std::vector<Pathway> path_from_start;
        for (auto point: extra_way_points)
            path_from_start.push_back((this->*algorithm)(start, point));

        std::vector<Pathway> path_between_points;
        for (auto point: extra_way_points)
            for (auto another_point: extra_way_points)
                if (point != another_point)
                    path_between_points.push_back((this->*algorithm)(point, another_point));

        std::vector<Pathway> path_to_destination;
        for (auto point: extra_way_points)
            path_to_destination.push_back((this->*algorithm)(point, destination));

        std::unordered_map<Node *, size_t> point_index;
        for (size_t i = 0; i < extra_way_points.size(); ++i)
            point_index[extra_way_points[i]] = i;
        std::set<size_t> points_left;
        for (size_t i = 1; i < extra_way_points.size() + 1; ++i)
            points_left.insert(i);

        return combinatory_pathway(points_left, 0,
                                   path_from_start, path_between_points, path_to_destination);
    }
};

int main() {
    // 60.000860, 30.368114 - дом
    double home_latitute = 60.000860;
    double home_longitute = 30.368114;
    // 59.944168, 30.295489 - итмо биржа
    double work_latitute = 59.944168;
    double work_longitute = 30.295489;

    My_way my_way(home_latitute, home_longitute,
                  work_latitute, work_longitute,
                  "7hw/spb_graph.txt");
    std::cout << my_way.start->latitute << " " << my_way.start->longitute << std::endl;
    std::cout << my_way.destination->latitute << " " << my_way.destination->longitute << std::endl;

    // 60.004211, 30.299131 - сити молл
    my_way.add_extra_point(60.004211, 30.299131);
    std::cout << my_way.extra_way_points.back()->latitute << " " << my_way.extra_way_points.back()->longitute <<
            std::endl;

    // 59.908491, 30.513576 - кудрово
    my_way.add_extra_point(59.908491, 30.513576);
    std::cout << my_way.extra_way_points.back()->latitute << " " << my_way.extra_way_points.back()->longitute <<
            std::endl;

    // 60.021335, 30.309469 - психбольница 3
    // my_way.add_extra_point(60.021335, 30.309469);
    // std::cout << my_way.extra_way_points.back()->latitute << " " << my_way.extra_way_points.back()->longitute <<
    //         std::endl;

    auto [path, dist] = my_way.find_best_path(&My_way::a_star);

    std::cout << "Path: ";
    for (auto node: path)
        if (node == my_way.destination)
            std::cout << node->latitute << " " << node->longitute << std::endl;
        else
            std::cout << node->latitute << " " << node->longitute << " -> " << std::endl;
    std::cout << std::endl;

    std::cout << "Distance: " << dist << std::endl;
}

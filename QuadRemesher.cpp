#include <vector>
#include <map>
#include <cstdint>
#include <cmath>
#include <limits>
#include <fstream>
#include <iostream>
#include <memory>
#include <array>

#include <vector>
#include <array>
#include <unordered_map>
#include <unordered_set>
#include <cstdint>

struct coord
{
    double x, y, z;
    coord() : x(0), y(0), z(0) {}
    coord(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
    coord operator+(const coord &other) const { return {x + other.x, y + other.y, z + other.z}; }
    coord operator-(const coord &other) const { return {x - other.x, y - other.y, z - other.z}; }
    coord operator*(double scalar) const { return {x * scalar, y * scalar, z * scalar}; }
    coord operator/(double scalar) const { return {x / scalar, y / scalar, z / scalar}; }
    double length() const { return std::sqrt(x * x + y * y + z * z); }
    double magnitude() const { return std::sqrt(x * x + y * y + z * z); }
    coord normalized() const { return *this / length(); }
    double abs() const { return std::abs(x) + std::abs(y) + std::abs(z); }

    double dot(const coord &other) const { return x * other.x + y * other.y + z * other.z; }
    coord cross(const coord &other) const { return {y * other.z - z * other.y, z * other.x - x * other.z, x * other.y - y * other.x}; }
};

// Forward declarations
class Quad;
class QuadNode;
class QuadMesh;

// Definition of the classes

/// @brief  A quad element
/// @details A quad element is defined by four QuadNode object pointers
class Quad
{
private:
    std::array<QuadNode *, 4> nodes;
    uint32_t id;
    double weight;
    // coord strain_vector;
public:
    explicit Quad(const std::array<QuadNode *, 4> &nodes, uint32_t in_id) : nodes(nodes), id(in_id), weight(0.0) {}

    // Setters
    void set_weight(double w) { weight = w; }

    // Getters
    std::array<QuadNode *, 2> get_adjacent_nodes(const QuadNode *p) const;
    const std::array<QuadNode *, 4>& get_nodes() const { return nodes; }

    uint32_t get_id() const { return id; }    
    double   get_weight() const { return weight; }
    double   get_area() const;
    coord    get_normal() const;
    coord    get_center() const 
    { return (nodes[0]->get_position() + nodes[1]->get_position() + nodes[2]->get_position() + nodes[3]->get_position()) / 4.0; }
};

/// @brief  A node in the mesh
/// @details A node is defined by its position and a list of quads that share this node
class QuadNode
{
private:
    coord position;
    std::vector<Quad *> belonging_quads;
    uint32_t id;
    double strain_factor;
    coord tension;

public:
    explicit QuadNode(const coord &p, uint32_t in_id) : position(p), id(in_id) , tension({0,0,0}) {}

    const coord& get_position() const { return position; }
    void set_position(const coord &p) { position = p; }

    void add_belonging_quad(Quad *quad) { belonging_quads.push_back(quad); }
    const std::vector<Quad *> &get_belonging_quads() const { return belonging_quads; }

    bool is_singularity() const { return get_connected_neighbour_nodes().size() % 2 == 1; }

    std::vector<QuadNode *> get_connected_neighbour_nodes() const;

    uint32_t get_id() const { return id; }
    const coord &get_strain() const { return tension * strain_factor; }
    const coord &get_tension_vector() const { return tension; }

    void set_strain_factor(const double factor) { strain_factor = factor; }
    void set_tension_vector(const coord &t) { tension = t; }
};


/// @brief  A mesh of nodes and quads
/// @details A mesh is defined by a list of nodes and a list of quads
class QuadMesh
{
private:
    std::vector<QuadNode> nodes;
    std::vector<Quad> quads;

public:
    QuadMesh(const std::vector<coord> &node_positions, const std::vector<std::array<uint32_t, 4>> &quad_faces);
    const std::vector<QuadNode> &get_nodes() const { return nodes; }
    const std::vector<Quad> &get_quads() const { return quads; }

    void print_info();
    void relax();

private:
    void generate_connectivity();
};

///********************************************************************************************************************///
///* Implementations                                                                                                  *///
///********************************************************************************************************************///

inline std::array<QuadNode *, 2> Quad::get_adjacent_nodes(const QuadNode* p) const
{
    if (!p)
    {
        return {nullptr, nullptr};
    }

    constexpr uint32_t adjacent_node_map[4][2] = {{1, 3}, {0, 2}, {1, 3}, {0, 2}};
    for (size_t i = 0; i < 4; ++i)
    {
        if (nodes[i] == p)
        {
            return {nodes[adjacent_node_map[i][0]], nodes[adjacent_node_map[i][1]]};
        }
    }
    return {nullptr, nullptr};
}

inline coord Quad::get_normal() const
{
    coord v1 = nodes[1]->get_position() - nodes[0]->get_position();
    coord v2 = nodes[3]->get_position() - nodes[0]->get_position();
    return v1.cross(v2).normalized();
}

inline double Quad::get_area() const
{
    coord v1 = nodes[1]->get_position() - nodes[0]->get_position();
    coord v2 = nodes[3]->get_position() - nodes[0]->get_position();
    coord v3 = nodes[2]->get_position() - nodes[0]->get_position();

    double area1 = 0.5 * v1.cross(v2).magnitude();
    double area2 = 0.5 * v3.cross(v2).magnitude();

    return area1 + area2;
}

inline std::vector<QuadNode *> QuadNode::get_connected_neighbour_nodes() const
{
    std::vector<QuadNode *> neighbours;
    std::unordered_set<QuadNode *> unique_neighbours;

    for (const Quad *quad : belonging_quads)
    {
        auto adjacent = quad->get_adjacent_nodes(this);
        for (QuadNode *adj_node : adjacent)
        {
            if (adj_node && unique_neighbours.insert(adj_node).second)
            {
                neighbours.push_back(adj_node);
            }
        }
    }
    return neighbours;
}

inline QuadMesh::QuadMesh(const std::vector<coord> &node_positions, const std::vector<std::array<uint32_t, 4>> &quad_faces)
{
    nodes.reserve(node_positions.size());
    uint32_t node_id = 0;
    for (const auto &pos : node_positions)
    {
        nodes.emplace_back(pos, node_id++);
    }

    quads.reserve(quad_faces.size());
    uint32_t quad_id = 0;
    for (const auto &face : quad_faces)
    {
        std::array<QuadNode *, 4> quad_nodes = {&nodes[face[0]], &nodes[face[1]], &nodes[face[2]], &nodes[face[3]]};
        quads.emplace_back(quad_nodes, quad_id++);
    }
    generate_connectivity();
}

inline void QuadMesh::generate_connectivity()
{
    for (Quad &quad : quads)
    {
        for (QuadNode *p : quad.get_nodes())
        {
            p->add_belonging_quad(&quad);
        }
    }
}

inline void QuadMesh::relax()
{
    for (QuadNode& node : nodes)
    {
        if(node.get_strain().magnitude() > 0.0)
        {
            std::cout << "Node " << node.get_id() << " has strain: " << node.get_strain().x << " " << node.get_strain().y << " " << node.get_strain().z << std::endl;
        }

        // Propagate the strain to the connected nodes & distribute it
        std::vector<QuadNode *> connected_nodes = node.get_connected_neighbour_nodes();

        coord new_center;
        
        for (QuadNode* p : connected_nodes)
        {
           new_center = new_center +  p->get_position();
        }

        new_center = new_center / connected_nodes.size();

        new_center = (new_center + node.get_position())*(0.5);

    }
}

inline void QuadMesh::print_info()
{
    for (const QuadNode& node : nodes)
    {
        std::cout << "Node position: " << node.get_position().x << " " << node.get_position().y << " " << node.get_position().z << std::endl;

        std::cout << "Belonging quads: ";
        for (Quad *p : node.get_belonging_quads())
        {
            std::cout << p->get_id() << " | ";
        }
        
        std::cout << std::endl;

        std::cout << "Connected nodes: ";

        for (QuadNode* p : node.get_connected_neighbour_nodes())
        {
            std::cout << "| " << p->get_id() << " | ";
        }
        std::cout << std::endl;
    }

    for (const Quad& quad : quads)
    {
        std::cout << "Quad nodes: ";
        for (QuadNode* p : quad.get_nodes())
        {
            std::cout << p->get_position().x << " " << p->get_position().y << " " << p->get_position().z << " | ";
        }
        std::cout << std::endl;
    }
}


int main()
{
    // Create a two quads sharing a common edge
    std::vector<coord> node_positions = {
        {0, 0, 0},
        {1, 0, 0},
        {1, 1, 0},
        {0, 1, 0},
        {1, 0, 1},
        {1, 1, 1},
        {0, 1, 1},
        {0, 0, 1},
    };

    std::vector<std::array<uint32_t, 4>> quad_faces = {
        {0, 1, 2, 3},
        {1, 4, 5, 2},
        {3, 2, 5, 6},
        {0, 3, 6, 7},
    };

    QuadMesh mesh(node_positions, quad_faces);
    mesh.print_info();
}

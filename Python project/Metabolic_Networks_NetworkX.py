import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx

class MetabolicNetwork ():
    
    def __init__(self, network_type = "metabolite-reaction", split_rev = False):
        self.G = nx.DiGraph()
        self.net_type = network_type
        self.node_types = {}
        if network_type == "metabolite-reaction":
            self.node_types["metabolite"] = []
            self.node_types["reaction"] = []
        self.split_rev = split_rev
    
    # Adds nodes to to the graph and the respective type (reaction or metabolite)
    def add_vertex_type(self, v, nodetype):
        self.G.add_node(v)
        self.node_types[nodetype].append(v)
    
    # Adds interactions between nodes 'o' and 'd'      
    def add_edge(self, o, d):
        self.G.add_edge(o, d)
    
    # Removes interactions between nodes 'u' and 'v'
    def remove_edge(self, u, v):
        nx.remove_edge(self.G, u, v)
    
    # Returns the type of node that exists in the graph
    def get_nodes_type(self, node_type):
        if node_type in self.node_types:
            return self.node_types[node_type]
        else: return None
        
    # Returns all the nodes in the graph
    def get_nodes(self):
        return self.G.nodes()

    # Returns all the interactions between all nodes
    def get_edges(self):
        return self.G.edges()
    
    # Returns the sucessors of node 'v'. Biologically, it returns the metabolites/reactions that depend on the node 'v'. 
    def get_successors(self, v):
        return self.G.sucessors(v)
    
    # Returns the predecessors of node 'v'. Biologically, it returns the metabolites/reactions that create node 'v'
    def get_predecessors(self, v):
        return self.G.predecessors(v)
    
    #Returns the out degree of a node, that is, the amount of reactions that consume that node
    def out_degree(self, v):
        return self.G.out_degree(v)
    
    # Returns the in degree of a node, that is, the amount of reactions that create the node
    def in_degree(self, v):
        return self.G.in_degree(v)
      
    # Returns the degree of a node 'v' (in and out).
    def degree(self, v):
        return self.G.degree(v)
    
    # Returns the degrees of all nodes in the graph
    def all_degrees(self):
        degrees = [(node,val) for (node, val) in self.G.degree()]
        return degrees
    
    # Returns the total average of degrees of all nodes in the graph
    def mean_degree (self):
        degrees = self.all_degrees()
        dg_int = [x[1] for x in degrees] 

        return sum(dg_int)/(len(dg_int))
    
    # Returns the shortest path between nodes 's' and 'd'
    def shortest_path(self, s, d):
        return nx.shortest_path(self.G,s,d)
    
    # Returns the distance of the shortest path between 's' and 'd'
    def distance(self, s, d):
        return nx.shortest_path_length(self.G,s,d)
    
    # Returns the average of the shortest paths of the graph
    def mean_shortest_path(self):
        return nx.average_shortest_path_length(self.G)
    
    # Verifies the existence of a path between two nodes 's' and 'd'
    def has_path(self,s,d):
        return nx.has_path(self.G,s,d)
    
    # Returns the average conectivity degree of a node, that is, returns the average degree of the directly linked nodes of a certain node with degree k
    def average_degree_connectivity(self):
        return nx.k_nearest_neighbors(self.G)
    
    # Returns the clustering coeficients of node 'v'
    def clustering_coef(self,v):
        return nx.clustering(self.G.to_undirected(),v)
    
    # Returns the clustering coeficients of all nodes in the graph
    def all_clustering_coefs(self):
        ccs = {}
        for k in list(self.G.nodes):
            ccs[k] = self.clustering_coef(k)
        
        return ccs
    
    # Returns the average of the clustering coeficients of all nodes in the graph
    def mean_clustering_coef(self):
        return nx.average_clustering(self.G.to_undirected())
    
    # Returns the hamiltonian path.
    def hamilton_path(self):
        return nx.hamilton_path(self.G.to_undirected())
    
    # Returns the edges of the eulerian circuit
    def eulerian_circuit(self):
        if nx.is_eulerian(self.G) == True:
            return list(nx.eulerian_circuit(self.G))
        else:
            raise Exception("Graph is not Eulerian!")  

    # Prints the graph with the shortest path between two different nodes 'source' and 'dest'
    def print_graph_sp(self,source,dest):

        sp = (m.shortest_path(source,dest))
        list_of_sp = []
        for i in range(len(sp)-1):
            list_of_sp.append((sp[i],sp[i+1]))
        
        red_edges = list_of_sp

        black_edges = [edge for edge in self.G.edges() if edge not in red_edges]
        
        pos = nx.spring_layout(self.G,7)
        
        metabolite_nodes = self.get_nodes_type("metabolite")
        reaction_nodes = self.get_nodes_type("reaction")
        
        nx.draw_networkx_nodes(self.G, pos, nodelist=metabolite_nodes, 
                               node_color = 'g', node_size = 500)
        
        nx.draw_networkx_nodes(self.G, pos, nodelist=reaction_nodes,
                               node_color = 'y', node_size = 500)
        
        nx.draw_networkx_labels(self.G, pos)
        nx.draw_networkx_edges(self.G, pos, edgelist=red_edges, edge_color='r', arrows=True)
        nx.draw_networkx_edges(self.G, pos, edgelist=black_edges, edge_color='black', arrows=True)
        
        plt.axis('off')
        plt.show()
    
    # Prints the resulting graph, distinguishing the reactions from the metabolites with different colours
    def print_graph(self):

        pos = nx.spring_layout(self.G,4)
        
        metabolite_nodes = self.get_nodes_type("metabolite")
        reaction_nodes = self.get_nodes_type("reaction")
        
        nx.draw_networkx_nodes(self.G, pos, nodelist=metabolite_nodes, 
                               node_color = 'g', node_size = 500)
        
        nx.draw_networkx_nodes(self.G, pos, nodelist=reaction_nodes, 
                               node_color = 'y', node_size = 500)
        
        nx.draw_networkx_labels(self.G, pos)
        nx.draw_networkx_edges(self.G, pos, edge_color='black', arrows=True)
        
        plt.axis('off')
        plt.show()

    
    # Read file
    def read_from_file(self,filename):
        
        f = open(filename,'r')
        finfo = f.read()
        f.close()

        term = str(filename).split('.')[1]  
        
        gmr = MetabolicNetwork("metabolite-reaction")
        
        if term == 'txt':
           
            reactions = finfo.split('\n')
            for line in reactions:
                tokens = line.split(':')
                reac_id = tokens[0].strip()
                gmr.add_vertex_type(reac_id, "reaction")
                rline = tokens[1]             
                if "<=>" in rline:
                    left, right = rline.split("<=>")
                    mets_left = left.split("+")
                    for met in mets_left:
                        met_id = met.strip()
                        if met_id not in gmr.G:
                            gmr.add_vertex_type(met_id, "metabolite")
                        if self.split_rev:
                            gmr.add_vertex_type(reac_id+"_b", "reaction")
                            gmr.add_edge(met_id, reac_id)
                            gmr.add_edge(reac_id+"_b", met_id)
                        else:
                            gmr.add_edge(met_id, reac_id)
                            gmr.add_edge(reac_id, met_id)
                    mets_right = right.split("+")
                    for met in mets_right:
                        met_id = met.strip()
                        if met_id not in gmr.G:
                            gmr.add_vertex_type(met_id, "metabolite")
                        if self.split_rev:
                            gmr.add_edge(met_id, reac_id+"_b")
                            gmr.add_edge(reac_id, met_id)
                        else:
                            gmr.add_edge(met_id, reac_id)
                            gmr.add_edge(reac_id, met_id)
                elif "=>" in line:
                    left, right = rline.split("=>")
                    mets_left = left.split("+")
                    for met in mets_left:
                        met_id = met.strip()
                        if met_id not in gmr.G:
                            gmr.add_vertex_type(met_id, "metabolite")
                        gmr.add_edge(met_id, reac_id)
                    mets_right = right.split("+")
                    for met in mets_right:
                        met_id = met.strip()
                        if met_id not in gmr.G:
                            gmr.add_vertex_type(met_id, "metabolite")
                        gmr.add_edge(reac_id, met_id)
    
        self.G = gmr.G
        self.node_types = gmr.node_types

if __name__ == '__main__':
    x = MetabolicNetwork("metabolite-reaction")
    
    x.read_from_file("example-net.txt")
    print(x.mean_shortest_path())
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
SEGMENTAL_EDGE = 'S'
ADJACENCY_EDGE = 'A'
SUPPLEMENTARY_EDGE_ID = 100000

import networkx as nx
from collections import Counter
import copy

def get_opposite(v):
    if v % 2 == 0:
        return v + 1
    return v - 1

def get_canonical_cycle(cycle):
    
    size = len(cycle)
    doubled_cycle = cycle[:-1] + cycle[:-1]
    reversed_doubled_cycle = list(reversed(doubled_cycle))

    return min([doubled_cycle[i:i + size] for i in range(size - 1)] + [reversed_doubled_cycle[i:i + size + 1] for i in range(size - 1)])


class Graph:

    def __init__(self):
        self.G = nx.Graph()
        
    def is_valid(self):
        
        return all([ self.get_copy_number_excess(v) >= 0 for v in self.G.nodes ])
        
    def is_linearly_decomposable(self):
         
        return self.is_valid() and nx.is_connected(self.G) and any([ self.get_copy_number_excess(v) > 0 for v in self.G.nodes ])
           
    def is_cyclically_decomposable(self):
    
        return self.is_valid() and nx.is_connected(self.G) and all([ self.get_copy_number_excess(v) == 0 for v in self.G.nodes ])
        
    def add_edge(self, v1, v2, edge_type, multiplicity=1):
        if not self.G.has_edge(v1, v2):
            self.G.add_edge(v1, v2)
            self.G[v1][v2][SEGMENTAL_EDGE] = 0
            self.G[v1][v2][ADJACENCY_EDGE] = 0
        self.G[v1][v2][edge_type] += multiplicity
            
        assert self.G[v1][v2][SEGMENTAL_EDGE] >= 0 
        
        if self.G[v1][v2][SEGMENTAL_EDGE] > 0:
            assert min(v1, v2) + 1 == max(v1, v2) and min(v1, v2) % 2 == 0
            
        assert self.G[v1][v2][ADJACENCY_EDGE] >= 0
        assert self.G[v1][v2][SEGMENTAL_EDGE] + self.G[v1][v2][SEGMENTAL_EDGE] >= 0   
                
        if  self.G[v1][v2][SEGMENTAL_EDGE] + self.G[v1][v2][ADJACENCY_EDGE] == 0:
            self.G.remove_edge(v1, v2)
                
    def add_segmental_edge(self, v, multiplicity=1):
        self.add_edge(v // 2 * 2, v // 2 * 2 + 1, SEGMENTAL_EDGE, multiplicity)
        
    def add_adjacency_edge(self, v1, v2, multiplicity=1):
        self.add_edge(v1, v2, ADJACENCY_EDGE, multiplicity)
        
    def remove_segmental_edge(self, v, multiplicity=1):
        self.add_segmental_edge(v, -multiplicity)

    def remove_adjacency_edge(self, v1, v2, multiplicity=1):
        self.add_adjacency_edge(v1, v2, -multiplicity)
        
    def get_copy_number_excess(self, v):
        adjacencies_multiplicity = sum(self.G[v][v2][ADJACENCY_EDGE] for v2 in self.G[v])
        
        if self.G.has_edge(v, v):
            adjacencies_multiplicity += self.G[v][v][ADJACENCY_EDGE]
        
        return self.get_segmental_edge_multiplicity(v) - adjacencies_multiplicity

    def is_telomere(self, v):
        return self.get_copy_number_excess(v) > 0
    
    def get_telomeres(self):
        return [ v for v in self.G.nodes if self.is_telomere(v) ]
    
    def get_overall_copy_number_excess(self):
        return sum([self.get_copy_number_excess(v) for v in self.G.nodes])

    def get_adjacency_neighbours(self, v):
        
        if len(self.G[v]) == 0: 
            return []
        return sorted([ v1 for v1 in self.G[v] if self.G[v][v1][ADJACENCY_EDGE] > 0])
    
    def get_segmental_edges(self):
        return [ sorted([v1, v2]) for (v1, v2, s) in self.G.edges(data=SEGMENTAL_EDGE) if s != 0 ] 
    
    def get_segmental_edge_multiplicity(self, v):
        
        s1 = v // 2 * 2
        s2 = v // 2 * 2 + 1
        
        return self.G[s1][s2][SEGMENTAL_EDGE] if self.G.has_edge(s1, s2) else 0
    
    def get_next_good_neighbour(self, v, last):

        adjacency_neighbours = self.get_adjacency_neighbours(v)
        adjacency_neighbours_next = [ v1 for v1 in  adjacency_neighbours if last is None or  v1 > last ][:2]
        
        if len(adjacency_neighbours_next) == 0:
            return None
        
        if v ==  SUPPLEMENTARY_EDGE_ID + 1 and last is not None:
            return None
        
        if last is None:
            return adjacency_neighbours[0]
        
        v2 = adjacency_neighbours_next[0]
            
        self.remove_adjacency_edge(v, v2)

        is_connected = nx.has_path(self.G, SUPPLEMENTARY_EDGE_ID, v)

        self.add_adjacency_edge(v, v2)            

        if is_connected:
            return v2
        
        if len(adjacency_neighbours_next) > 1:
            return adjacency_neighbours_next[1]
                                          
    def get_multiplicities_sum(self):
        return sum(self.G[v1][v2][SEGMENTAL_EDGE] + self.G[v1][v2][ADJACENCY_EDGE] for v1, v2 in self.G.edges)
    
    def is_empty(self):
        return len(self.G.edges) == 0                                                                                            
        
    def build_connectivity_certificate(self):

        telomeres = self.get_telomeres()
        k = self.get_overall_copy_number_excess() // 2

        connectivity_certificate = copy.deepcopy(self)
        connectivity_certificate.add_segmental_edge(SUPPLEMENTARY_EDGE_ID, 2 * k)

        if k > 1:
            connectivity_certificate.add_adjacency_edge(SUPPLEMENTARY_EDGE_ID, SUPPLEMENTARY_EDGE_ID, k - 1)

        for v in telomeres:
            x = self.get_copy_number_excess(v)
            connectivity_certificate.add_adjacency_edge(SUPPLEMENTARY_EDGE_ID + 1, v, x)

        return connectivity_certificate

    def enumerate_moeds(self):

        def split_connectivity_certificate_trail(trail):

            splitted_trail = [[]]

            for i in range(2, len(trail)):
                if trail[i - 1] == SUPPLEMENTARY_EDGE_ID and  trail[i - 2] == SUPPLEMENTARY_EDGE_ID + 1:
                    splitted_trail.append([])

                if trail[i] not in (SUPPLEMENTARY_EDGE_ID, SUPPLEMENTARY_EDGE_ID + 1):
                    splitted_trail[-1].append(trail[i])

            return splitted_trail
        
        connectivity_certificate = self.build_connectivity_certificate()
        
        trail = [ SUPPLEMENTARY_EDGE_ID, SUPPLEMENTARY_EDGE_ID + 1 ]
        connectivity_certificate.remove_segmental_edge(SUPPLEMENTARY_EDGE_ID)

        last_adjacency_vertex = None

        while True:

            next_good_neighbour = connectivity_certificate.get_next_good_neighbour(trail[-1], last_adjacency_vertex)

            if next_good_neighbour is not None:
                connectivity_certificate.remove_adjacency_edge(trail[-1], next_good_neighbour)
                connectivity_certificate.remove_segmental_edge(next_good_neighbour)
                trail.append(next_good_neighbour)
                trail.append(next_good_neighbour + (1 if next_good_neighbour % 2 == 0  else -1))

                last_adjacency_vertex = None

                if connectivity_certificate.is_empty():
                     yield split_connectivity_certificate_trail(trail)

            else:
                if len(trail) == 2:   
                    break

                connectivity_certificate.add_segmental_edge(trail[-1])
                connectivity_certificate.add_adjacency_edge(trail[-3], trail[-2])

                last_adjacency_vertex = trail[-2]
                trail = trail[:-2]
    
    def __deepcopy__(self, memodict={}):
        new_instance = Graph()
        new_instance.G = self.G.copy()
        return new_instance

    def __str__(self):
        
        text = ''
        
        for v1, v2  in self.G.edges:
            text += f"{v1} {v2}: {self.G[v1][v2]}\n"
        
        return text
    
    def get_reduced_graph(self, trails):
        
        kg = copy.deepcopy(self)
        
        for trail in trails:
            for i, (v1, v2) in enumerate(zip(trail, trail[1:])):
                if i % 2 == 0:
                    kg.remove_segmental_edge(v1)
                else:
                    kg.remove_adjacency_edge(v1, v2)   

        return kg
    
    def generate_next_trails(self):
    
        def generate_trails(kg, prefix, telomeres):

            v = prefix[-1]

            if v in telomeres:
                yield prefix

            for u in kg.get_adjacency_neighbours(v):   
                kg.remove_adjacency_edge(v, u)
                kg.remove_segmental_edge(u)
                yield from generate_trails(kg, prefix + [u, get_opposite(u)],  telomeres) 
                kg.add_adjacency_edge(v, u)
                kg.add_segmental_edge(u)              

        telomere = sorted(self.get_telomeres())[0]
        self.remove_segmental_edge(telomere)
        telomeres = self.get_telomeres()
        yield from generate_trails(self, [telomere , get_opposite(telomere)], telomeres)
        self.add_segmental_edge(telomere)

    def generate_next_cycles(self):    

        def generate_cycles(kg, prefix):

            v = prefix[-1]
            start_edge = prefix[0]
            
            for u in kg.get_adjacency_neighbours(v):   
                
                if start_edge == u:
                    yield prefix + [u]

                if not (u == start_edge and self.get_segmental_edge_multiplicity(start_edge) == 0):

                    kg.remove_adjacency_edge(v, u)
                    kg.remove_segmental_edge(u)
               
                    yield from generate_cycles(kg, prefix + [u, get_opposite(u)]) 
                    kg.add_adjacency_edge(v, u)
                    kg.add_segmental_edge(u)              
                    
        start_vertex = sorted(self.get_segmental_edges())[0][0]
        self.remove_segmental_edge(start_vertex)
        yield from generate_cycles(self, [start_vertex, get_opposite(start_vertex)])
        self.add_segmental_edge(start_vertex)
        
    def enumerate_walks(self, generator, stop_condition):
        
        prefixes =  [[ walk ] for walk in generator(self)() ]

        while prefixes:
            prefix = prefixes.pop()
            
            kg = self.get_reduced_graph(prefix)
      
            if stop_condition(kg):
                yield prefix
            else:
                for walk in generator(kg)():
                    prefixes.append(prefix + [walk])
        
    def enumerate_cycles(self):
        
        yield from self.enumerate_walks(lambda kg: kg.generate_next_cycles, lambda kg: kg.is_empty()) 
    
    def enumerate_trails(self):
        
        yield from self.enumerate_walks(lambda kg: kg.generate_next_trails, lambda kg: len(kg.get_telomeres()) == 0) 
                    
    def enumerate_naively(self):
        
        for trails in self.enumerate_trails():
            kg = self.get_reduced_graph(trails)
            if kg.is_empty():
                yield trails
            else:
                for cycles in kg.enumerate_cycles():
                  
                    yield trails +  sorted([get_canonical_cycle(cycle) for cycle in cycles])
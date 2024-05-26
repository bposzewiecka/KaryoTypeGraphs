import random

random.seed(666)

def get_random_karyotype_graph_type1(size):

    g = Graph()

    cum_size = 0
    
    vertices = []
    
    for size in sizes:
        for i in range(cum_size, cum_size + size):
            g.add_segmental_edge(2 * i)
            if i != cum_size + size - 1:
                g.add_adjacency_edge(2 * i + 1, 2 * i + 2)
                
            if (i - cum_size ) % 2 == 1:
 

        cum_size += size
        
    random.shuffle(vertices)

    for i in range(len(vertices)//2):
        g.add_adjacency_edge(vertices[2 * i], vertices[2 * i + 1]) 
        
    assert g.is_valid()

    return g


def get_random_karyotype_graph_type2(sizes):

    g = Graph()

    cum_size = 0
    
    vertices = []
    
    edges = []
    
    for size in sizes:
        for i in range(cum_size, cum_size + size):
            g.add_segmental_edge(2 * i)
            if i != cum_size + size - 1:
                g.add_adjacency_edge(2 * i + 1, 2 * i + 2)                    
                
        cum_size += size
    
    additional_segmental_edges = [ random.choice(range(size)) for i in range(2 * size) ]
    
    for i in additional_segmental_edges:
        g.add_segmental_edge(2 * i)
        vertices.append(2 * i)
        vertices.append(2 * i + 1)
       
    random.shuffle(vertices)

    for i in range(len(vertices)//2):
        g.add_adjacency_edge(vertices[2 * i], vertices[2 * i + 1]) 

    assert g.is_valid()
    
    return g

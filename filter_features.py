from scipy import stats
import sys

class FNode:
    def __init__(self, id, values, ks):
        self.feature_id = id
        self.values = values
        self.ks = ks
    def __eq__(self, other):
        return self.feature_id == other.feature_id

    def __str__(self):
        return str(self.feature_id)

class FEdge:
    def __init__(self, ids, r):
        self.ids = ids
        self.r = r

    def __str__(self):
        return str(self.r)

class Subgraphs:
    def __init__(self, fl):
        self.subgraphs = []
        self.features_list = fl

    def add(self, c):
        self.subgraphs.append(c)

    def __str__(self):
        res = ''

        res += str(len(self.subgraphs)) + ' subgraphs'

        for sg in self.subgraphs:
            nodes = sg[:]
            nodes.sort(key=lambda x: int(str(x)), reverse=False)
            res += '\n' + str(len(sg)) + ' nodes: ' + str([str(n) for n in nodes])
            res += '\n' + str(len(sg)) + ' nodes: ' + str([self.features_list[int(str(n))] for n in nodes])


            res += '\n'

        return res

def parse_csv_to_feature_list(filename, feature_start, feature_end):
    # id, data1, data2, ..., class
    graph = []
    features = []
    print(filename)
    with open(filename) as fp:
        line = fp.readline()
        features_list = line.rstrip("\n").split(',')[feature_start:feature_end]


        line = fp.readline()
        feature_list = line.rstrip("\n").split(',')[feature_start:feature_end]
        for feature in feature_list:
            features.append([float(feature)])

        while line:
            line = fp.readline()
            feature_list = line.rstrip("\n").split(',')[feature_start:feature_end]
            for i in range(len(feature_list)):
                features[i].append(float(feature_list[i]))

    return features, features_list

def feature_lists_to_nodes(feature_lists):
    nodes = []

    for i in range(len(feature_lists)):
        feature_values = feature_lists[i]
        ks = stats.kstest(feature_values, 'kstwobign', alternative = 'two-sided')
        nodes.append(FNode(i, feature_values, ks.statistic))

    return nodes

def filter_nodes(nodes):

    filtered = [x for x in nodes if x.ks > 0.1]
    filtered.sort(key=lambda x: x.ks, reverse=False)

    if (len(filtered) > 19):
        filtered = filtered[:19]

    return filtered

def nodes_to_edges(nodes):

    edges = []

    for i in range(len(nodes)):
        for j in range(i+1, len(nodes)):
            r = stats.pearsonr(nodes[i].values, nodes[j].values)[0]
            edges.append(FEdge((nodes[i].feature_id, nodes[j].feature_id), r))

    return edges

def filter_edges(edges):

    edges.sort(key=lambda x: abs(x.r), reverse=True)
    print([str(c) for c in edges])

    filtered = edges

    filtered = filtered[int(0.099 * len(filtered)):]
    filtered = [x for x in filtered if abs(x.r) <= 0.9]
    return filtered

def find_max_id(nodes):
    max_id = 0
    for n in nodes:
        max_id = max(n.feature_id, max_id)
    return max_id

def edges_to_adjacency_matrix(edges, size):
    adj = [0] * size
    adj = [adj.copy() for i in range(size)]
    for edge in edges:
        i = edge.ids[0]
        j = edge.ids[1]
        adj[i][j] = edge
        adj[j][i] = edge

    return adj

def nodes_to_map(nodes):
    map = {}
    for n in nodes:
        map[n.feature_id] = n
    return map

def bk_mcl(cliques, edges, potential_clique=[], remaining_nodes=[], skip_nodes=[], depth=0):

    # To understand the flow better, uncomment this:
    # print (' ' * depth), 'potential_clique:', potential_clique, 'remaining_nodes:', remaining_nodes, 'skip_nodes:', skip_nodes

    if len(remaining_nodes) == 0 and len(skip_nodes) == 0:
        cliques.add(potential_clique)
        return

    for node in remaining_nodes[:]:

        # Try adding the node to the current potential_clique to see if we can make it work.
        new_potential_clique = potential_clique + [node]
        new_remaining_nodes = [n for n in remaining_nodes if edges[node.feature_id][n.feature_id] != 0]
        new_skip_list = [n for n in skip_nodes if edges[node.feature_id][n.feature_id] != 0]
        bk_mcl(cliques, edges, new_potential_clique, new_remaining_nodes, new_skip_list, depth + 1)

        # We're done considering this node.  If there was a way to form a clique with it, we
        # already discovered its maximal clique in the recursive call above.  So, go ahead
        # and remove it from the list of remaining nodes and add it to the skip list.
        remaining_nodes.remove(node)
        skip_nodes.append(node)

def evaluate_sg(sg, adj_mat):
    edge_sum = 0
    edge_count = 0
    for i in range(len(sg)):
        for j in range(i+1, len(sg)):
            if (adj_mat[sg[i].feature_id][sg[j].feature_id] != 0):
                edge = adj_mat[sg[i].feature_id][sg[j].feature_id]
                edge_sum += abs(edge.r)
                edge_count += 1
    edge_sum = edge_sum/edge_count

    node_sum = 0
    for n in sg:
        node_sum += n.ks
    node_sum = node_sum/len(sg)

    return 1 - edge_sum + node_sum

def find_best_sg(subgraphs, adj_mat):
    max = evaluate_sg(subgraphs[0], adj_mat)
    best_sg = subgraphs[0]
    for sg in subgraphs[1:]:
        score = evaluate_sg(sg, adj_mat)
        if (score > max):
            max = score
            best_sg = sg
    return best_sg

def main():
    filename = sys.argv[1]
    fstart = sys.argv[2]
    fend = sys.argv[3]
    features, features_list = parse_csv_to_feature_list(filename, int(fstart), int(fend))
    nodes = feature_lists_to_nodes(features)
    filtered_nodes = filter_nodes(nodes)
    edges = nodes_to_edges(filtered_nodes)
    filtered_edges = filter_edges(edges)
    node_map = nodes_to_map(filtered_nodes)
    max_id = find_max_id(filtered_nodes)
    adj_mat = edges_to_adjacency_matrix(filtered_edges, max_id + 1)

    mcl = Subgraphs(features_list)
    bk_mcl(mcl, adj_mat, remaining_nodes=filtered_nodes)
    #
    print(mcl)

    best = find_best_sg(mcl.subgraphs, adj_mat)
    print(len(best))

    best.sort(key=lambda x: int(str(x)), reverse=False)

    for b in best:
        print(b, features_list[int(str(b))])


main()

import re

def default_separation(node_a, node_b):
    if node_a.parent == node_b.parent:
        return 1
    return 2

def equal_separation(node_a, node_b):
    return 1

class Tree(object):
    def __init__(self, ID = None, children = list(), x = 0, y = 0, 
                 label = None, length = 0, depth = None, parent = None,
                 cumulative_length = 0):
        self.ID = ID
        self.children = children
        self.x = x
        self.y = y
        self.label = label
        self.depth = depth
        self.parent = parent
        self.length = length
        self.cumulative_length = cumulative_length
    def __repr__(self):
        return f'<TreeNode ID={self.ID} depth={self.depth} length={self.length}>'
    
    @classmethod
    def from_newick(cls, newick_string):
        tokens = re.split('\s*(;|\(|\)|,|:)\s*', newick_string)
        ID = 0
        tree = cls(ID = ID, length = 0, cumulative_length = 0)
        ancestors = list()
        for i,token in enumerate(tokens):
            if token == '(':
                ID += 1
                subtree = cls(ID = ID)
                tree.children = [subtree]
                ancestors.append(tree)
                tree = subtree
            elif token == ',':
                ID += 1
                subtree = cls(ID = ID)
                ancestors[-1].children.append(subtree)
                tree = subtree
            elif token == ')':
                tree = ancestors.pop()
            else:
                previous_token = tokens[i - 1]
                if previous_token in ('(',')',','):
                    tree.name = token
                elif previous_token == ':':
                    tree.length = float(token)
        tree.depth = 0
        queue = [tree]
        while queue:
            node = queue.pop(0)
            for child in node.children:
                child.parent = node
                child.depth = node.depth + 1
                child.cumulative_length = node.cumulative_length + abs(child.length)
            queue += node.children
        
        return tree
    
    def to_newick(self, branch_lengths = True):
        if self.children:
            subtree_string = ','.join([
                c.to_newick(branch_lengths = branch_lengths) for c in self.children
            ])
            newick = f'({subtree_string}){self.name}'
        else:
            newick = self.name
        
        if branch_lengths and self.ID != 0:
            length = self.length
            if length == 0:
                length = int(0)
            newick += f':{length}'
        
        if self.ID == 0:
            newick += ';'
        
        return newick
    
    @property
    def nodes(self):
        return list(self.breadth_first())
    
    @property
    def leafs(self):
        return [n for n in self.nodes if not n.children]
    
    @property
    def links(self):
        _links = []
        for node in self.nodes:
            if node.children:
                for child in node.children:
                    _links.append((node, child))
        return _links
    
    @property
    def left_outer_leaf(self):
        node = self
        while node.children:
            node = node.children[0]
        return node
    
    @property
    def right_outer_leaf(self):
        node = self
        while node.children:
            node = node.children[-1]
        return node
    
    def breadth_first(self):
        queue = [self]
        while queue:
            node = queue.pop(0)
            queue += node.children
            yield node
        
    def depth_first(self, post_order = True):
        if not post_order:
            yield self
        for child in self.children:
            yield from child.depth_first(post_order = post_order)
        if post_order:
            yield self
    
    def layout(self, separation = default_separation, d_x = 1, d_y = 1, 
               ltr = True):
        previous_node = None
        y = 0
        for node in self.depth_first(post_order = True):
            if node.children:
                node.y = sum([c.y for c in node.children]) / len(node.children)
                if ltr:
                    node.x = 1 + max([c.x for c in node.children]) 
                else:
                    node.x = min([c.x for c in node.children]) - 1
            else:
                if previous_node:
                    y += separation(node, previous_node)
                    node.y = y
                else:
                    node.y = 0
                node.x = 0
                previous_node = node

        for node in self.depth_first(post_order = True):
            node.x = (self.x - node.x) * d_x
            node.y = (node.y - self.y) * d_y
        return self
